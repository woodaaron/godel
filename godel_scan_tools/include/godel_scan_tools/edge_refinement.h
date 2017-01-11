/*
    TODO:
    - Add more comments describing what each section of code is doing.
    - Move function implementations to separate cpp file.
    - Check ENSENSO scan density and predict the amount of neighbors.
    - SPEED THIS UP!!!
*/

#ifndef EDGE_REFINEMENT_H
#define EDGE_REFINEMENT_H

#include <math.h>
#include <pcl/common/common.h>
#include <pcl_ros/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d.h>
#include <pcl/surface/gp3.h>
#include <pcl/features/boundary.h>
#include <pcl/surface/concave_hull.h>
#include <pcl/kdtree/kdtree_flann.h>

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
  std::size_t current_pose_index_;
  std::size_t num_poses_;
  pcl::visualization::PCLVisualizer *viewer_;

  EigenPoseMatrix boundary_poses_;
  PointCloudVector boundary_pose_neighbor_;
  PointCloudVector refined_boundary_pose_neighbor_;
  PointCloudVector neighbor_boundary_points_;
  PointVector new_pose_points_;

  DebugDisplayData(const std::size_t current_pose_index, const std::size_t num_poses, 
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
    debug_display_ = false;
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

    for (std::size_t i = 0; i < boundary_poses.size(); i++)
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
        for (std::size_t j = 0; j < pointIdxRadiusSearch.size(); j++)
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

    for (std::size_t i = 0; i < boundary_poses.size(); i++)
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
        for (std::size_t j = 0; j < pointIdxNKNSearch.size(); j++)
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
    for (std::size_t i = 0; i < boundary_poses.size(); i++)
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
      for (std::size_t j = 0; j < boundary_pose_neighbor[i].size(); j++)
      {
        float x = boundary_pose_neighbor[i].points[j].x;
        float y = boundary_pose_neighbor[i].points[j].y;
        float z = boundary_pose_neighbor[i].points[j].z;

        float plane = calculatePlane(a, b, c, x, y, z, dot_product);
        deviations.push_back(std::abs(plane));        
      }

      float allowed_error = calculateAllowedDeviation(deviations);

      pcl::PointCloud<pcl::PointXYZ> temp_cloud;

      for (std::size_t j = 0; j < boundary_pose_neighbor[i].size(); j++)
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

    for (std::size_t i = 0; i < refined_cloud.size(); i++)
    {
      boundaries.clear();
      normals.clear();
      computeNormals(refined_cloud[i].makeShared(), normals);

      pcl::BoundaryEstimation<pcl::PointXYZ, pcl::Normal, pcl::Boundary> boundary_estimation;
      pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());

      boundary_estimation.setInputCloud(refined_cloud[i].makeShared());
      boundary_estimation.setInputNormals(normals.makeShared());
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

    for (std::size_t i = 0; i < refined_points_cloud.size(); i++)
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
    for (std::size_t i = 0; i < 4; i++)
    {
      for (std::size_t j = 0; j < 4; j++)
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
    for (std::size_t i = 0; i < original_boundary_poses.size(); i++)
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

    for (std::size_t i = 0; i < boundary_poses.size(); i++)
    {
      kdtree.setInputCloud(extracted_boundary_points[i].makeShared());

      searchpoint.x = boundary_poses[i](0, 3);
      searchpoint.y = boundary_poses[i](1, 3);
      searchpoint.z = boundary_poses[i](2, 3);

      pcl::PointXYZ temp_point;

      if (kdtree.nearestKSearch(searchpoint, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
      {
        for (std::size_t j = 0; j < pointIdxNKNSearch.size(); j++)
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

    for (std::size_t i = 0; i < boundary_poses.size(); i++)
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

    for (std::size_t i = 0; i < boundary_poses.size(); i++)
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

  void setDebugDisplay(bool debug_display)
  {
    if (debug_display) { debug_display_ = true; }
  }

  static void 
  keyboardEventOccurred(const pcl::visualization::KeyboardEvent &event,
                        void* debug_display_data_void)
  {
    DebugDisplayData *debug_display_data = static_cast<DebugDisplayData *> (debug_display_data_void);

    if (event.getKeySym() == "Right" && event.keyDown())
    {
      if (debug_display_data->current_pose_index_ >= 0 || debug_display_data->current_pose_index_ <= debug_display_data->num_poses_)
      {
        if (debug_display_data->current_pose_index_ == debug_display_data->num_poses_) { debug_display_data->current_pose_index_ = 0; }
        else { debug_display_data->current_pose_index_++; }
      }
    }

    if (event.getKeySym() == "Left" && event.keyDown())
    {
      if (debug_display_data->current_pose_index_ >= 0 || debug_display_data->current_pose_index_ <= debug_display_data->num_poses_)
      {
        if (debug_display_data->current_pose_index_ == 0) { debug_display_data->current_pose_index_ = debug_display_data->num_poses_; }
        else { debug_display_data->current_pose_index_--; }
      }
    }

    std::string display_text;
    display_text = "Current Pose: " + std::to_string(debug_display_data->current_pose_index_);
    debug_display_data->viewer_->updateText(display_text, 0, 15, "current pose");
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
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color_3(debug_display_data->neighbor_boundary_points_[debug_display_data->current_pose_index_].makeShared(), 255, 255, 0);
    debug_display_data->viewer_->addPointCloud<pcl::PointXYZ> (debug_display_data->neighbor_boundary_points_[debug_display_data->current_pose_index_].makeShared(), single_color_3, "Boundary Points");
    debug_display_data->viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1.5, "Boundary Points");
    // New Point
    debug_display_data->viewer_->addSphere(debug_display_data->new_pose_points_[debug_display_data->current_pose_index_], 2.5, 0.0, 1.0, 0.0, "new point");    
  }

  void
  debugDisplay(const EigenPoseMatrix &boundary_poses,
               const PointCloudVector &boundary_pose_neighbor,
               const PointCloudVector &refined_boundary_pose_neighbor,
               const PointCloudVector &neighbor_boundary_points,
               const PointVector &new_pose_points,
               const EigenPoseMatrix &refined_poses)
  {
    
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(visual_cloud_);
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer ("Debug Viewer"));
    viewer->setBackgroundColor(0, 0, 0);
    viewer->addCoordinateSystem(1.0);

    viewer->initCameraParameters();
    viewer->addPointCloud<pcl::PointXYZRGB> (visual_cloud_, "input cloud");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "input cloud");
    viewer->addText("Current Pose: ", 0, 15, "current pose");

    for (std::size_t i = 0; i < (boundary_poses.size() - 1); i++)
    {
      assert(boundary_poses.size() == refined_poses.size()); // May need to remove this in the future?
      std::string original_name = "original_line_" + std::to_string(i);
      std::string refined_name = "refined_line_" + std::to_string(i);
      pcl::PointXYZ original_p1(boundary_poses[i](0,3), boundary_poses[i](1,3), boundary_poses[i](2,3));
      pcl::PointXYZ original_p2(boundary_poses[i+1](0,3), boundary_poses[i+1](1,3), boundary_poses[i+1](2,3));
      pcl::PointXYZ refined_p1(refined_poses[i](0,3), refined_poses[i](1,3), refined_poses[i](2,3));
      pcl::PointXYZ refined_p2(refined_poses[i+1](0,3), refined_poses[i+1](1,3), refined_poses[i+1](2,3));

      viewer->addLine<pcl::PointXYZ>(original_p1, original_p2, original_name);
      viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, original_name);
      viewer->addLine<pcl::PointXYZ>(refined_p1, refined_p2, refined_name);
      viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, refined_name);
    }

    DebugDisplayData debug_display_data(current_pose_index_, num_poses_, viewer.get(), 
                                        boundary_poses, boundary_pose_neighbor, refined_boundary_pose_neighbor,
                                        neighbor_boundary_points, new_pose_points);

    viewer->registerKeyboardCallback(keyboardEventOccurred, static_cast<void *>(&debug_display_data));

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

  /*
      The returned index is the first index of the jump. So if the index is 100, the large jump is
      between poses 100 and 101.
  */

  static float
  distanceBetweenTwoPoints(const PointVector &point_vector,
                           const int index_1,
                           const int index_2)
  {
    float diff_x = point_vector[index_1].x - point_vector[index_2].x;
    float diff_y = point_vector[index_1].y - point_vector[index_2].y;
    float diff_z = point_vector[index_1].z - point_vector[index_2].z;

    float magnitude = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z);

    return magnitude;
  }

  static int
  calculateNumberOfPointsToInsert(const float &distance_between_points,
                                  const float &standard_deviation)
  {
    int number_of_points = round(distance_between_points / standard_deviation);
    return number_of_points - 2; // -2 because of the original two points.
  }

  static void
  calculateOutliersInNewPosePoints(const PointVector &neighbor_new_pose_points,
                                   std::map<int, int> &outlier_index)
  {
    std::vector<float> difference_between_poses;
    for (std::size_t i = 1; i < neighbor_new_pose_points.size(); i++)
    {
      float magnitude = distanceBetweenTwoPoints(neighbor_new_pose_points, i, i-1);
     
      if (magnitude != 0.0)
      {
        difference_between_poses.push_back(magnitude);
      }
    }

    float sum = std::accumulate(difference_between_poses.begin(), difference_between_poses.end(), 0.0);
    float mean = sum / difference_between_poses.size();
    float standard_deviation = calculateAllowedDeviation(difference_between_poses);

    std::cout << "Mean: " << mean << std::endl;
    std::cout << "Deviation: " << standard_deviation << std::endl;
    std::cout << "Max Deviation: " << maxValueOfVector(difference_between_poses) << std::endl;

    for (std::size_t i = 1; i < neighbor_new_pose_points.size(); i++)
    {
      float magnitude = distanceBetweenTwoPoints(neighbor_new_pose_points, i, i-1);

      if ((magnitude-3*standard_deviation) >= mean)
      {
        std::cout << "Pose: " << (i-1) << " - " << i << " : " << magnitude << std::endl;
        outlier_index[i-1] = calculateNumberOfPointsToInsert(magnitude, standard_deviation);
      } 
    }

    // Debug check
    for (std::map<int, int>::const_iterator it = outlier_index.begin(); it != outlier_index.end(); it++)
    {
      std::cout << "Pose: " << it->first << " requires " << it->second << " points." << std::endl;
    }
  }

  static PointVector
  calculateClosestBoundaryPointToNextPose(const PointCloudVector &boundary_points, 
                                          const PointVector &neighbor_new_pose_points, 
                                          const int &index,
                                          const int &num_poses_required)
  {
    PointVector additional_points;
    int K = 1;
    int pose_index;
    int closest_pose_index;

    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(boundary_points[index].makeShared());
    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);

    pcl::PointXYZ searchpoint = neighbor_new_pose_points[index];
    if (kdtree.nearestKSearch(searchpoint, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
    {
      {
        pose_index = pointIdxNKNSearch[0];
        std::cout << "Pose Index: " << pose_index << std::endl;
      }
    }

    if (kdtree.nearestKSearch(neighbor_new_pose_points[index+1], K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
    {
      closest_pose_index = pointIdxNKNSearch[0];
      std::cout << "Closest Pose Index: " << closest_pose_index << std::endl;
    }
    
    if (pose_index == 0 && closest_pose_index > 0)
    {
      // This means the shortest way to the other pose is to move to the next pose in the boundary.
      if ((closest_pose_index - pose_index) < boundary_points[index].width / 2)
      {
        std::cout << "The closest way to get to closest_pose_index is to add" << std::endl;
        if (closest_pose_index > num_poses_required)
        {
          int skip_point = (int)floor(closest_pose_index / num_poses_required);
          for (std::size_t i = 0; i < closest_pose_index; i += skip_point)
          {
            additional_points.push_back(boundary_points[index].points[i]);
          }
        }
        else
        {
          for (std::size_t i = 0; i < closest_pose_index; i++)
          {
            additional_points.push_back(boundary_points[index].points[i]);
          }
        }
      }
      // This means the shortest way to the other pose is to move backwards in the boundary.
      else
      {
        std::cout << "The closest way to get to the closest_pose_index is to subtract" << std::endl;
        if ((closest_pose_index - ((boundary_points[index].width / 2)-1)) > num_poses_required)
        {
          int skip_point = (int)floor((closest_pose_index - ((boundary_points[index].width / 2) - 1)) / num_poses_required);
          for (std::size_t i = 0; i < (closest_pose_index - ((boundary_points[index].width / 2) - 1)); i+= skip_point)
          {
            if (i == 0) { additional_points.push_back(boundary_points[index].points[i]); }
            else
            {
              int point_index = boundary_points[index].width - i - skip_point;
              additional_points.push_back(boundary_points[index].points[point_index]);
            }
          }
        }
        else
        {
          for (std::size_t i = 0; i < (closest_pose_index - ((boundary_points[index].width / 2) - 1)); i++)
          {
            if (i == 0) { additional_points.push_back(boundary_points[index].points[i]); }
            else
            {
              int point_index = boundary_points[index].width - i;
              additional_points.push_back(boundary_points[index].points[point_index]);
            }
          }
        }
      }
    }
    return additional_points;
  }

  static void
  iCantThinkOfANameForThisYet(const PointCloudVector &boundary_points, 
                              const PointVector &neighbor_new_pose_points, 
                              const std::map<int, int> &outlier_index,
                              std::map<int, PointVector> &additional_poses)
  {
    for (std::map<int, int>::const_iterator it = outlier_index.begin(); it != outlier_index.end(); it++)
    {
      int index = it->first;
      int num_poses_required = it->second;
      std::cout << "Index #: " << index << std::endl;
      additional_poses[it->first] = calculateClosestBoundaryPointToNextPose(boundary_points, 
                                      neighbor_new_pose_points, index, num_poses_required);
    }

    // Debug check
    for (std::map<int, PointVector>::const_iterator it = additional_poses.begin(); it != additional_poses.end(); it++)
    {
      std::cout << "Pose: " << it->first << " requires the following additional points:" << std::endl;
      for (std::size_t i = 0; i < it->second.size(); i++)
      {
        std::cout << it->second[i] << std::endl;
      }
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
    num_poses_ = boundary_poses.size() - 1;
    current_pose_index_ = num_poses_ / 2; // Sets the debug visualization so it starts near the top.

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

    // 5) Find any outliers in the new pose points (new points that jump too far).
    std::map<int, int> outlier_index; // Pose Number, Number of Points to Add
    calculateOutliersInNewPosePoints(neighbor_new_pose_points, outlier_index);

    // 6) Determines the boundary points that follow the shortest distance between the poses of the two outliers.
    std::map<int, PointVector> additional_poses;
    iCantThinkOfANameForThisYet(neighbor_boundary_points, neighbor_new_pose_points, outlier_index, additional_poses);

    // 7) Move original boundary pose point to new point while keeping same orientation
    // movePoseToNewPoint(boundary_poses, radius_new_pose_points, refined_poses);
    movePoseToNewPoint(boundary_poses, neighbor_new_pose_points, refined_poses);

    if (debug_display_)
    {
      debugDisplay(boundary_poses, boundary_pose_neighbor, refined_boundary_pose_neighbor, 
                   neighbor_boundary_points, neighbor_new_pose_points, refined_poses);
    }

    #if 0
    for (std::size_t i = 0; i < boundary_poses.size(); i++)
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

  static float maxValueOfVector(std::vector<float> &vec)
  {
    float max_value = 0;
    for (std::size_t i = 0; i < vec.size(); i++)
    {
      if (vec[i] > max_value)
      {
        max_value = vec[i];
      }
    }
    return max_value;
  }

private:
  pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud_;
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_;

  bool debug_display_;
  float point_density_;
  double radius_;
  double sradius_;
  int edge_direction_;

  float boundary_search_radius_;
  float search_radius_;
  int number_of_neighbors_;
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr visual_cloud_;

  std::size_t current_pose_index_;
  std::size_t num_poses_;
};
} // namespace godel_scan_tools
#endif // EDGE_REFINEMENT_H