/*
    TODO:
    - Add more comments describing what each section of code is doing.
    [x] Calculate tolerance of a local plane by checking the standard deviation of the deviation 
      of the nearby points to determine if it is an outlier.
    - Adjust boundary estimator constant.
    [?] Set boundary threshold value. (Set the angle threshold for a point to be a normal)
    [1/2] Define a variable for boundary threshold and comment regarding the physical meaning.
    - Create a debugging display for updating one pose.
      - Show original point cloud in one color, the points within a radius in another color,
        those that fit on the plane in another, the boundary points in another.
        - Do this display step by step to make sure it is doing exactly what we expect at each step.
*/

#ifndef EDGE_REFINEMENT_H
#define EDGE_REFINEMENT_H

#include <functional>
#include <pcl/common/common.h>
#include <pcl_ros/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d.h>
#include <pcl/surface/gp3.h>
#include <pcl/features/boundary.h>
#include <pcl/surface/concave_hull.h>
#include <pcl/kdtree/kdtree_flann.h>

// Yea this is bad but i'm too tired to figure out how to implement this.
size_t CURRENT_POSE_INDEX = 0;
size_t NUM_POSES;

namespace godel_scan_tools
{

typedef std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> PointCloudVector;
typedef std::vector<pcl::PointCloud<pcl::Boundary>, Eigen::aligned_allocator<pcl::Boundary>> PointCloudBoundaryVector;
typedef std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> EigenPoseMatrix;
typedef std::vector<pcl::PointXYZ> PointVector;
typedef Eigen::Matrix<float, 1, 3> NormalVector;
typedef Eigen::Matrix<float, 1, 3> PoseOrigin;

struct DebugDisplayData
{
  size_t current_pose_index_;
  size_t num_poses_;
  pcl::visualization::PCLVisualizer *viewer_;

  EigenPoseMatrix boundary_poses_;
  PointCloudVector boundary_pose_neighbor_;
  PointCloudVector refined_boundary_pose_neighbor_;
  PointCloudVector neighbor_boundary_points_;
  PointVector new_pose_points_;

  DebugDisplayData(const size_t current_pose_index, const size_t num_poses, 
                   pcl::visualization::PCLVisualizer *viewer,
                   const EigenPoseMatrix boundary_poses, 
                   const PointCloudVector boundary_pose_neighbor, 
                   const PointCloudVector refined_boundary_pose_neighbor, 
                   const PointCloudVector neighbor_boundary_points,
                   const PointVector new_pose_points)
  {
    current_pose_index_ = current_pose_index;
    num_poses_ = num_poses;
    viewer_ = viewer;

    boundary_poses_ = boundary_poses;
    boundary_pose_neighbor_ = boundary_pose_neighbor;
    refined_boundary_pose_neighbor_ = refined_boundary_pose_neighbor;
    neighbor_boundary_points_ = neighbor_boundary_points;
    new_pose_points_ = new_pose_points;
  }
};

class EdgeRefinement
{
public:
  /*
  *   @brief constructor 
  *   @param cloud input cloud from which you plan to refine a boundary
  */
  EdgeRefinement(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud):
                tree_(new pcl::search::KdTree<pcl::PointXYZ>()), 
                point_density_(0), edge_direction_(0)
  {
    tree_->setInputCloud(cloud);
    input_cloud_= pcl::PointCloud<pcl::PointXYZ>::Ptr(cloud);
    visual_cloud_ = pcl::PointCloud<pcl::PointXYZRGB>::Ptr(new pcl::PointCloud<pcl::PointXYZRGB>);
    current_pose_index_ = 0;
    getPointDensity();
  }

#if 1
  float 
  getPointDensity(void)
  {
    if(point_density_ > 0) 
      return(point_density_);

    int n = input_cloud_->points.size();
    int K = 100;

    pcl::PointXYZ pt = input_cloud_->points[n/1];
    std::vector<int> pt_indices(K);
    std::vector<float> pt_sq_distances(K);

    int num_found;
    if((num_found = tree_->nearestKSearch(pt, K, pt_indices, pt_sq_distances))>0)
    {
      double maxd = 0;
      int maxi = 0;
      for(int i = 0; i < K; i++) // note, there should be one point with zero distance
      {
        if(maxd < pt_sq_distances[i]) 
        {
          maxd = pt_sq_distances[i];
          maxi = i;
        }
      }

      //printf("maxd = %lf\n",maxd);

      double r = sqrt(maxd);
      double v = 4/3*3.14*r*r*r; /* volume containing K points Kpts/V */

      point_density_ = K/v; // k pts per vol
      radius_ = cbrt(3/(4*3.13*point_density_))/150;

      //printf("calculated radius_ = %f\n",radius_);
      sradius_ = radius_;
    }

    else
    {
      printf("Could not find %d points near center of input_cloud_\n",K);
      point_density_ = 0.0;
    }
    return(point_density_);
  }
#endif

  static void
  nearestNeighborRadiusSearch(const pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud,
                              const EigenPoseMatrix &boundary_poses,
                              PointCloudVector &boundary_pose_neighbors,
                              const float search_radius)
  {
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(input_cloud);

    for (size_t i = 0; i < boundary_poses.size(); i++)
    {
      std::vector<int> pointIdxRadiusSearch;
      std::vector<float> pointRadiusSquaredDistance;

      pcl::PointXYZ searchpoint;
      searchpoint.x = boundary_poses[i](0, 3);
      searchpoint.y = boundary_poses[i](1, 3);
      searchpoint.z = boundary_poses[i](2, 3);

      pcl::PointCloud<pcl::PointXYZ> temp_cloud;
      if (kdtree.radiusSearch(searchpoint, search_radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
      {
        for (size_t j = 0; j < pointIdxRadiusSearch.size(); j++)
        {
          temp_cloud.push_back(input_cloud->points[pointIdxRadiusSearch[j]]);
        }
      }

      boundary_pose_neighbors.push_back(temp_cloud);
    }
  }

  static void 
  nearestNNeighborSearch(const pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud,
                         const EigenPoseMatrix &boundary_poses,
                         PointCloudVector &boundary_pose_neighbor,
                         const int n)
  {
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(input_cloud);

    for (size_t i = 0; i < boundary_poses.size(); i++)
    {
      std::vector<int> pointIdxNKNSearch(n);
      std::vector<float> pointNKNSquaredDistance(n);

      pcl::PointXYZ searchpoint;
      searchpoint.x = boundary_poses[i](0, 3);
      searchpoint.y = boundary_poses[i](1, 3);
      searchpoint.z = boundary_poses[i](2, 3);

      pcl::PointCloud<pcl::PointXYZ> temp_cloud;
      if (kdtree.nearestKSearch(searchpoint, n, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
      {
        for (size_t j = 0; j < pointIdxNKNSearch.size(); j++)
        {
          temp_cloud.push_back(input_cloud->points[pointIdxNKNSearch[j]]);
        }
      }

      boundary_pose_neighbor.push_back(temp_cloud);
    }
  }

  /*
      a*x + b*y + c*z = d
  */
  static float
  calculatePlane(const float &x, const float &y, const float &z, 
                 const float &a, const float &b, const float &c, 
                 const float &d)
  {
    return ((a*x + b*y + c*z) - d);
  }

  /* 
      To check if a point is in a plane:
      normal.x * x + normal.y * y + normal.z * z = origin dot normal
      normal.x * x + normal.y * y + normal.z * z - (origin dot normal) = 0 
  */
  static void 
  refineNeighborPoints(const EigenPoseMatrix &boundary_poses,
                       const PointCloudVector &boundary_pose_neighbor,
                       PointCloudVector &refined_boundary_pose_neighbor)
  {
    for (size_t i = 0; i < boundary_poses.size(); i++)
    { 
      NormalVector normal;
      PoseOrigin pose_origin;

      std::vector<float> deviations;
      deviations.reserve(boundary_pose_neighbor[i].size());      

      normal(0, 0) = boundary_poses[i](0, 2);
      normal(0, 1) = boundary_poses[i](1, 2);
      normal(0, 2) = boundary_poses[i](2, 2);

      pose_origin(0, 0) = boundary_poses[i](0, 3);
      pose_origin(0, 1) = boundary_poses[i](1, 3);
      pose_origin(0, 2) = boundary_poses[i](2, 3);

      float a = normal(0, 0);
      float b = normal(0, 1);
      float c = normal(0, 2);

      float dot_product = pose_origin.dot(normal);

      // Calcualtes the deviation of the nearby points from the plane.
      for (size_t j = 0; j < boundary_pose_neighbor[i].size(); j++)
      {
        float x = boundary_pose_neighbor[i].points[j].x;
        float y = boundary_pose_neighbor[i].points[j].y;
        float z = boundary_pose_neighbor[i].points[j].z;

        float plane = calculatePlane(a, b, c, x, y, z, dot_product);
        deviations.push_back(std::abs(plane));        
      }

      float allowed_error = calculateAllowedDeviation(deviations);

      pcl::PointCloud<pcl::PointXYZ> temp_cloud;

      for (size_t j = 0; j < boundary_pose_neighbor[i].size(); j++)
      {
        float x = boundary_pose_neighbor[i].points[j].x;
        float y = boundary_pose_neighbor[i].points[j].y;
        float z = boundary_pose_neighbor[i].points[j].z;

        float plane = calculatePlane(a, b, c, x, y, z, dot_product);

        if (plane <= allowed_error && plane >= -allowed_error)
        {
          temp_cloud.push_back(boundary_pose_neighbor[i].points[j]);
        }
      }

      refined_boundary_pose_neighbor.push_back(temp_cloud);
    }
  }

  // Calcualtes standard deviation of the deviations.
  // http://stackoverflow.com/questions/7616511/calculate-mean-and-standard-deviation-from-a-vector-of-samples-in-c-using-boos
  static float 
  calculateAllowedDeviation(const std::vector<float> &deviations)
  {
    float allowed_deviation;

    float sum = std::accumulate(deviations.begin(), deviations.end(), 0.0);
    float mean = sum / deviations.size();

    std::vector<float> diff(deviations.size());
    std::transform(deviations.begin(), deviations.end(), diff.begin(), [mean](double x) { return x - mean; });
    float sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    float stdev = std::sqrt(sq_sum / deviations.size());

    allowed_deviation = stdev;

    return allowed_deviation;
  }

  static void
  computeNormals(const pcl::PointCloud<pcl::PointXYZ>::Ptr &input_cloud, pcl::PointCloud<pcl::Normal> &normals)
  {
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normal_estimation;
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());

    normal_estimation.setInputCloud(input_cloud);
    normal_estimation.setSearchMethod(tree);
    normal_estimation.setKSearch(5);
    normal_estimation.compute(normals);
  }

  static void
  computeBoundaryForRefinedCloud(const PointCloudVector &refined_cloud,
                                 PointCloudBoundaryVector &refined_boundary,
                                 const float boundary_search_radius)
  {
    pcl::PointCloud<pcl::Boundary> boundaries;
    pcl::PointCloud<pcl::Normal> normals;

    for (size_t i = 0; i < refined_cloud.size(); i++)
    {
      pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud (new pcl::PointCloud<pcl::PointXYZ>);
      for (const auto &pt : refined_cloud[i].points)
      {
        input_cloud->push_back(pt);
      }

      boundaries.clear();
      normals.clear();
      computeNormals(input_cloud, normals);
      pcl::PointCloud<pcl::Normal>::Ptr normal_ptr (new pcl::PointCloud<pcl::Normal>);

      for (const auto &n : normals.points)
      {
        normal_ptr->push_back(n);
      }

      pcl::BoundaryEstimation<pcl::PointXYZ, pcl::Normal, pcl::Boundary> boundary_estimation;
      pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());

      boundary_estimation.setInputCloud(input_cloud);
      boundary_estimation.setInputNormals(normal_ptr);
      boundary_estimation.setRadiusSearch(boundary_search_radius);
      boundary_estimation.setSearchMethod(tree);
      boundary_estimation.setAngleThreshold(90.0 * 3.14 / 180.0); // Defaults to PI/2 according to the documentation...
      boundary_estimation.compute(boundaries);

      refined_boundary.push_back(boundaries);
    }
  }

  static void
  extractBoundaryPointsFromPointCloud(const PointCloudVector &refined_points_cloud,
                                      const PointCloudBoundaryVector &boundary_cloud,
                                      PointCloudVector &boundary_points)
  {
    const float bad_point = std::numeric_limits<float>::quiet_NaN();

    pcl::PointCloud<pcl::PointXYZ> temp_cloud;

    for (size_t i = 0; i < refined_points_cloud.size(); i++)
    {
      temp_cloud.clear();
      int k = 0;
      for (const auto &pt : boundary_cloud[i].points)
      {
        if (pt.boundary_point)
        {
          temp_cloud.push_back(refined_points_cloud[i].points[k]);
        }

        k++;
      }

      boundary_points.push_back(temp_cloud);
    }     
  }

  static bool
  containsNaNs(Eigen::Matrix4f matrix)
  {
    for (size_t i = 0; i < 4; i++)
    {
      for (size_t j = 0; j < 4; j++)
      {
        if (std::isnan(matrix(i, j))) { return true; }
        else { return false; }
      }
    }
  }

  static void 
  removeNaNFromPoseTrajectory(const EigenPoseMatrix &original_boundary_poses,
                              EigenPoseMatrix &boundary_poses_no_nan)
  {
    for (size_t i = 0; i < original_boundary_poses.size(); i++)
    {
      if (!containsNaNs(original_boundary_poses[i]))
      {
        boundary_poses_no_nan.push_back(original_boundary_poses[i]);
      }
    }
  }

  static void
  calculateClosestPointInBoundaryToPose(const EigenPoseMatrix &boundary_poses,
                                        const PointCloudVector &extracted_boundary_points,
                                        PointVector &new_pose_points)
  {
    int K = 1;
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    pcl::PointXYZ searchpoint;

    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);

    for (size_t i = 0; i < boundary_poses.size(); i++)
    {
      kdtree.setInputCloud(extracted_boundary_points[i].makeShared());

      searchpoint.x = boundary_poses[i](0, 3);
      searchpoint.y = boundary_poses[i](1, 3);
      searchpoint.z = boundary_poses[i](2, 3);

      pcl::PointXYZ temp_point;

      if (kdtree.nearestKSearch(searchpoint, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
      {
        for (size_t j = 0; j < pointIdxNKNSearch.size(); j++)
        {
          temp_point.x = extracted_boundary_points[i].points[pointIdxNKNSearch[j]].x;
          temp_point.y = extracted_boundary_points[i].points[pointIdxNKNSearch[j]].y;
          temp_point.z = extracted_boundary_points[i].points[pointIdxNKNSearch[j]].z;
        }

        new_pose_points.push_back(temp_point);
      }
    }
  }

  static void
  comparePoints(const EigenPoseMatrix &boundary_poses,
                const PointVector &radius_new_boundary_points,
                const PointVector &neighbor_new_boundary_points)
  {
    assert(boundary_poses.size() == radius_new_boundary_points.size());
    assert(boundary_poses.size() == neighbor_new_boundary_points.size());

    for (size_t i = 0; i < boundary_poses.size(); i++)
    {
      std::cout << "Old Boundary Point: " << "x: " << boundary_poses[i](0, 3) 
                << ", y: " << boundary_poses[i](1, 3) << ", z: " 
                << boundary_poses[i](2, 3) << std::endl;
      std::cout << "Radius New Boundary Point: " << "x: " 
                << radius_new_boundary_points[i].x << ", y: " 
                << radius_new_boundary_points[i].y 
                << ", z: " << radius_new_boundary_points[i].z << std::endl;
      std::cout << "Neighbor New Boundary Point: " << "x: " 
                << neighbor_new_boundary_points[i].x << ", y: " 
                << neighbor_new_boundary_points[i].y 
                << ", z: " << neighbor_new_boundary_points[i].z << std::endl << std::endl;                
    }
  }

  static void
  movePoseToNewPoint(const EigenPoseMatrix &boundary_poses,
                     const PointVector &new_boundary_points,
                     EigenPoseMatrix &refined_poses)
  {
    Eigen::Matrix4f temp_pose;

    for (size_t i = 0; i < boundary_poses.size(); i++)
    {
      temp_pose = boundary_poses[i];
      temp_pose(0, 3) = new_boundary_points[i].x;
      temp_pose(1, 3) = new_boundary_points[i].y;
      temp_pose(2, 3) = new_boundary_points[i].z;

      refined_poses.push_back(temp_pose);
    }
  }

  void 
  setBoundarySearchRadius(float search_radius)
  {
    boundary_search_radius_ = search_radius;
  }

  float
  getBoundarySearchRadius(void)
  {
    return boundary_search_radius_;
  }

  void 
  setNumberOfNeighbors(int number_of_neighbors)
  {
    number_of_neighbors_ = number_of_neighbors;
  }

  int 
  getNumberOfNeighbors(void)
  {
    return number_of_neighbors_;
  }

  void 
  setSearchRadius(float search_radius)
  {
    search_radius_ = search_radius;
  }

  float 
  getSearchRadius(void)
  {
    return search_radius_;
  }

  static void 
  keyboardEventOccurred(const pcl::visualization::KeyboardEvent &event,
                        void* debug_display_data_void)
                        //void* viewer_void)
  {
    //pcl::visualization::PCLVisualizer *viewer = static_cast<pcl::visualization::PCLVisualizer *> (viewer_void);
    DebugDisplayData *debug_display_data = static_cast<DebugDisplayData *> (debug_display_data_void);

    if (event.getKeySym() == "Right" && event.keyDown())
    {
      if (debug_display_data->current_pose_index_ >= 0 || debug_display_data->current_pose_index_ <= debug_display_data->num_poses_)
      {
        if (debug_display_data->current_pose_index_ == debug_display_data->num_poses_) { debug_display_data->current_pose_index_ = 0; }
        else { debug_display_data->current_pose_index_++; }
      }
      std::cout << debug_display_data->current_pose_index_ << std::endl;
    }

    if (event.getKeySym() == "Left" && event.keyDown())
    {
      if (debug_display_data->current_pose_index_ >= 0 || debug_display_data->current_pose_index_ <= debug_display_data->num_poses_)
      {
        if (debug_display_data->current_pose_index_ == 0) { debug_display_data->current_pose_index_ = debug_display_data->num_poses_; }
        else { debug_display_data->current_pose_index_--; }
      }
      std::cout << debug_display_data->current_pose_index_ << std::endl;
    }

    debug_display_data->viewer_->removeShape("pose point");
    debug_display_data->viewer_->removeShape("new point");
    debug_display_data->viewer_->removePointCloud("nearest N neighbors");
    debug_display_data->viewer_->removePointCloud("N neighbors in plane");
    debug_display_data->viewer_->removePointCloud("Boundary Points");

    pcl::PointXYZ pose_point;
    pose_point.x = debug_display_data->boundary_poses_[debug_display_data->current_pose_index_](0, 3);
    pose_point.y = debug_display_data->boundary_poses_[debug_display_data->current_pose_index_](1, 3);
    pose_point.z = debug_display_data->boundary_poses_[debug_display_data->current_pose_index_](2, 3);
    // Pose Point
    debug_display_data->viewer_->addSphere(pose_point, 2.5, 1.0, 0.0, 0.0, "pose point");
    // Points within certain radius or K neighbors
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color_1(debug_display_data->boundary_pose_neighbor_[debug_display_data->current_pose_index_].makeShared(), 0, 255, 0);
    debug_display_data->viewer_->addPointCloud<pcl::PointXYZ> (debug_display_data->boundary_pose_neighbor_[debug_display_data->current_pose_index_].makeShared(), single_color_1, "nearest N neighbors");
    debug_display_data->viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "nearest N neighbors");
    // Points within plane
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color_2(debug_display_data->refined_boundary_pose_neighbor_[debug_display_data->current_pose_index_].makeShared(), 0, 0, 255);
    debug_display_data->viewer_->addPointCloud<pcl::PointXYZ> (debug_display_data->refined_boundary_pose_neighbor_[debug_display_data->current_pose_index_].makeShared(), single_color_2, "N neighbors in plane");
    debug_display_data->viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "N neighbors in plane");
    // Boundary Points
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color_3(debug_display_data->neighbor_boundary_points_[debug_display_data->current_pose_index_].makeShared(), 255, 0, 0);
    debug_display_data->viewer_->addPointCloud<pcl::PointXYZ> (debug_display_data->neighbor_boundary_points_[debug_display_data->current_pose_index_].makeShared(), single_color_3, "Boundary Points");
    debug_display_data->viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "Boundary Points");
    // New Point
    debug_display_data->viewer_->addSphere(debug_display_data->new_pose_points_[debug_display_data->current_pose_index_], 2.5, 255.0/255.0, 69.0/255.0, 0.0, "new point");    
  }

  void
  debugDisplay(const EigenPoseMatrix &boundary_poses,
               const PointCloudVector &boundary_pose_neighbor,
               const PointCloudVector &refined_boundary_pose_neighbor,
               const PointCloudVector &neighbor_boundary_points,
               const PointVector &new_pose_points)
  {
    
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(visual_cloud_);
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer ("Debug Viewer"));
    viewer->setBackgroundColor(0, 0, 0);
    viewer->addCoordinateSystem(1.0);

    viewer->initCameraParameters();
    viewer->addPointCloud<pcl::PointXYZRGB> (visual_cloud_, "input cloud");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "input cloud");

    DebugDisplayData debug_display_data(current_pose_index_, num_poses_, viewer.get(), 
                                        boundary_poses, boundary_pose_neighbor, refined_boundary_pose_neighbor,
                                        neighbor_boundary_points, new_pose_points);

    // viewer->registerKeyboardCallback(keyboardEventOccurred, (void*)viewer.get());
    viewer->registerKeyboardCallback(keyboardEventOccurred, static_cast<void *>(&debug_display_data));

#if 0
    size_t temp_i = 100;
    //for (size_t i = 0; i < boundary_poses.size(); i++)
    {
      pcl::PointXYZ pose_point;
      pose_point.x = boundary_poses[temp_i](0, 3);
      pose_point.y = boundary_poses[temp_i](1, 3);
      pose_point.z = boundary_poses[temp_i](2, 3);
      // Pose Point
      viewer->addSphere(pose_point, 2.5, 1.0, 0.0, 0.0, "pose point");
      // Points within certain radius or K neighbors
      pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color_1(boundary_pose_neighbor[temp_i].makeShared(), 0, 255, 0);
      viewer->addPointCloud<pcl::PointXYZ> (boundary_pose_neighbor[temp_i].makeShared(), single_color_1, "nearest N neighbors");
      viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "nearest N neighbors");
      // Points within plane
      pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color_2(refined_boundary_pose_neighbor[temp_i].makeShared(), 0, 0, 255);
      viewer->addPointCloud<pcl::PointXYZ> (refined_boundary_pose_neighbor[temp_i].makeShared(), single_color_2, "N neighbors in plane");
      viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "N neighbors in plane");
      // Boundary Points
      pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color_3(neighbor_boundary_points[temp_i].makeShared(), 255, 0, 0);
      viewer->addPointCloud<pcl::PointXYZ> (neighbor_boundary_points[temp_i].makeShared(), single_color_3, "Boundary Points");
      viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "Boundary Points");
      // New Point
      viewer->addSphere(new_pose_points[temp_i], 2.5, 255.0/255.0, 69.0/255.0, 0.0, "new point");
    }
#endif

    while (!viewer->wasStopped())
    {
      viewer->spinOnce (100);
      boost::this_thread::sleep(boost::posix_time::microseconds (100000));
    }    
  }

  void
  setVisualCloud(const pcl::PointCloud<pcl::PointXYZRGB>::Ptr colored_cloud_ptr)
  {
    visual_cloud_->clear();

    for (const auto &pt : *colored_cloud_ptr)
    {
      visual_cloud_->push_back(pt);
    }
  }

  void 
  refineBoundary(const EigenPoseMatrix &original_boundary_poses, 
                 EigenPoseMatrix &refined_poses)
  {
    refined_poses.clear();

    // Remove NaNs from input boundary poses.
    EigenPoseMatrix boundary_poses;
    boundary_poses.reserve(original_boundary_poses.size());
    removeNaNFromPoseTrajectory(original_boundary_poses, boundary_poses);
    num_poses_ = boundary_poses.size();

    // 1) Find all points within R1 of each boundary pose.
    PointCloudVector boundary_pose_radius;
    PointCloudVector boundary_pose_neighbor;
    boundary_pose_radius.reserve(boundary_poses.size());
    boundary_pose_neighbor.reserve(boundary_poses.size());

    // nearestNeighborRadiusSearch(input_cloud_, boundary_poses, boundary_pose_radius, search_radius_);
    nearestNNeighborSearch(input_cloud_, boundary_poses, boundary_pose_neighbor, number_of_neighbors_);

    // 2) Narrow down the radius points at each pose to only lie on the x-y plane of the pose with some error.
    PointCloudVector refined_boundary_pose_radius;
    PointCloudVector refined_boundary_pose_neighbor;

    refined_boundary_pose_radius.reserve(boundary_poses.size());
    refined_boundary_pose_neighbor.reserve(boundary_poses.size());

    // refineNeighborPoints(boundary_poses, boundary_pose_radius, refined_boundary_pose_radius);
    refineNeighborPoints(boundary_poses, boundary_pose_neighbor, refined_boundary_pose_neighbor);

    // 3) Find all points that are boundaries.
    PointCloudBoundaryVector radius_boundary;
    PointCloudBoundaryVector neighbor_boundary;

    radius_boundary.reserve(boundary_poses.size());
    neighbor_boundary.reserve(boundary_poses.size());
    
    // computeBoundaryForRefinedCloud(refined_boundary_pose_radius, radius_boundary, boundary_search_radius_);
    computeBoundaryForRefinedCloud(refined_boundary_pose_neighbor, neighbor_boundary, boundary_search_radius_);

    PointCloudVector radius_boundary_points;
    PointCloudVector neighbor_boundary_points;

    radius_boundary_points.reserve(boundary_poses.size());
    neighbor_boundary_points.reserve(boundary_poses.size());

    // extractBoundaryPointsFromPointCloud(refined_boundary_pose_radius, radius_boundary, radius_boundary_points);
    extractBoundaryPointsFromPointCloud(refined_boundary_pose_neighbor, neighbor_boundary, neighbor_boundary_points);

    // 4) Find the boundary point that is closest to the original.
    PointVector radius_new_pose_points;
    PointVector neighbor_new_pose_points;

    radius_new_pose_points.reserve(boundary_poses.size());
    neighbor_new_pose_points.reserve(boundary_poses.size());

    // calculateClosestPointInBoundaryToPose(boundary_poses, radius_boundary_points, radius_new_pose_points);
    calculateClosestPointInBoundaryToPose(boundary_poses, neighbor_boundary_points, neighbor_new_pose_points);

    // 5) Move original boundary pose point to new point while keeping same orientation
    // movePoseToNewPoint(boundary_poses, radius_new_pose_points, refined_poses);
    movePoseToNewPoint(boundary_poses, neighbor_new_pose_points, refined_poses);

    debugDisplay(boundary_poses, boundary_pose_neighbor, refined_boundary_pose_neighbor, 
                 neighbor_boundary_points, neighbor_new_pose_points);

    #if 0
    for (size_t i = 0; i < boundary_poses.size(); i++)
    {
      std::cout << std::endl << "Boundary Pose Number: " << i << std::endl;

      std::cout << "(Radius Search) Number of Neighbors: " << boundary_pose_radius[i].width << std::endl;
      std::cout << "(Radius Search) Number of Refined Neighbors: " << refined_boundary_pose_radius[i].width << std::endl;
      std::cout << "(Radius Search) Number of Boundary Solutions: " << radius_boundary[i].width << std::endl;
      std::cout << "(Radius Search) Number of Boundary Points: " << radius_boundary_points[i].width << std::endl;

      std::cout << "(N-Neighbor Search) Number of Neighbors: " << boundary_pose_neighbor[i].width << std::endl;
      std::cout << "(N-Neighbor Search) Number of Refined Neighbors: " << refined_boundary_pose_neighbor[i].width << std::endl;  
      std::cout << "(N-Neighbor Search) Number of Boundary Solutions: " << neighbor_boundary[i].width << std::endl;  
      std::cout << "(N-Neighbor Search) Number of Boundary Points: " << neighbor_boundary_points[i].width << std::endl;

      std::cout << std::endl;        
    }
    #endif
  }

private:
  pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud_;
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_;

  float point_density_;
  double radius_;
  double sradius_;
  int edge_direction_;

  float boundary_search_radius_;
  float search_radius_;
  int number_of_neighbors_;
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr visual_cloud_;

  size_t current_pose_index_;
  size_t num_poses_;
};
} // namespace godel_scan_tools
#endif // EDGE_REFINEMENT_H