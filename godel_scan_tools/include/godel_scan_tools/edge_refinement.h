/*
    TODO:
    - Add more comments describing what each section of code is doing.
    - Move function implementations to separate cpp file.
    - Check ENSENSO scan density and predict the amount of neighbors.
    - SPEED THIS UP!!!
    - Add new points to refined_poses and possibly check if order is correct.
    - Add B-Spline Smoother:
      http://stackoverflow.com/questions/25379422/b-spline-curves
      http://kluge.in-chemnitz.de/opensource/spline/
      http://pointclouds.org/documentation/tutorials/bspline_fitting.php
*/

#ifndef EDGE_REFINEMENT_H
#define EDGE_REFINEMENT_H

#include <math.h>
#include <pcl/common/common.h>
#include <pcl_ros/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/surface/gp3.h>
#include <pcl/features/boundary.h>
#include <pcl/surface/concave_hull.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/visualization/pcl_visualizer.h>      

namespace godel_scan_tools
{
typedef std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> PointCloudVector;
typedef std::vector<pcl::PointCloud<pcl::Boundary>, Eigen::aligned_allocator<pcl::Boundary>> PointCloudBoundaryVector;
typedef std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> EigenPoseMatrix;
typedef std::vector<pcl::PointXYZ> PointVector;
typedef Eigen::Matrix<float, 1, 3> NormalVector;
typedef Eigen::Matrix<float, 1, 3> PoseOrigin;

/**
 * @brief      Structure containing the data for the debug display.
 */
struct DebugDisplayData
{
  /**
   * @brief      Constructor for DebugDisplayData.
   *
   * @param[in]  current_pose_index              The current pose index
   * @param[in]  num_poses                       The number poses
   * @param      viewer                          The viewer
   * @param[in]  boundary_poses                  The boundary poses
   * @param[in]  boundary_pose_neighbor          The boundary pose neighbor
   * @param[in]  refined_boundary_pose_neighbor  The refined boundary pose neighbor
   * @param[in]  neighbor_boundary_points        The neighbor boundary points
   * @param[in]  new_pose_points                 The new pose points
   * @param[in]  additional_poses                The additional poses
   */
  DebugDisplayData(const std::size_t current_pose_index, const std::size_t num_poses, 
                   pcl::visualization::PCLVisualizer *viewer,
                   const EigenPoseMatrix boundary_poses, 
                   const PointCloudVector boundary_pose_neighbor, 
                   const PointCloudVector refined_boundary_pose_neighbor, 
                   const PointCloudVector neighbor_boundary_points,
                   const PointVector new_pose_points,
                   const std::map<int, PointVector> additional_poses);

  bool rendered_additional_shapes_;
  std::size_t rendered_shape_count_;
  std::size_t current_pose_index_;
  std::size_t num_poses_;
  pcl::visualization::PCLVisualizer *viewer_;

  EigenPoseMatrix boundary_poses_;
  PointCloudVector boundary_pose_neighbor_;
  PointCloudVector refined_boundary_pose_neighbor_;
  PointCloudVector neighbor_boundary_points_;
  PointVector new_pose_points_;
  std::map<int, PointVector> additional_poses_;

};

/**
 * @brief      Class for edge refinement.
 */
class EdgeRefinement
{
public:
  /**
   * @brief      Constructor for EdgeRefinement class.
   *
   * @param[in]  cloud  Original Point cloud data that does not contain any NaNs.
   */
  EdgeRefinement(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud);

  /**
   * @brief      Constructor for EdgeRefinement class with initializer list.
   */
  EdgeRefinement():tree_(new pcl::search::KdTree<pcl::PointXYZ>()), point_density_(0), edge_direction_(0) {}

  /**
   * @brief      Main function of the class. Refines the boundary defined by the original boundary
   *             poses and returns the refined poses.
   *
   * @details    This algorithm removes any NaN's from the original boundary poses, then:
   *             1) Finds all points within N points of each boundary pose using a nearest neighbor search.
   *             2) Removes all points that do not lie on the x-y plane of that pose with some allowed deviation.
   *             3) Determines which of the remaining points are "boundary points".
   *             4) Extracts those boundary point indices from the original point cloud.
   *             5) Determines which boundary point is closest to the original pose.
   *             6) Generates refined poses by moving the original pose orientations to the closest boundary point.
   *             7) Calculates if any of the new pose points make a large jump.
   *             8) Generates additional points that follow along the path of the actual boundary for the large jump.
   *             8) Determines which indices require addional poses and adds them to the refined poses.
   *              
   * @param[in]  original_boundary_poses  The original boundary poses
   * @param      refined_poses            The refined poses
   */
  void refineBoundary(const EigenPoseMatrix &original_boundary_poses, EigenPoseMatrix &refined_poses);

  /**
   * @brief      Sets the search radius for the boundary search.
   *
   * @param[in]  search_radius  The search radius
   */
  void setBoundarySearchRadius(float search_radius) { boundary_search_radius_ = search_radius; }

  /**
   * @brief      Gets the search radius for the boundary search.
   *
   * @return     The boundary search radius.
   */
  float getBoundarySearchRadius(void) { return boundary_search_radius_; }

  /**
   * @brief      Sets the number of neighbors to find for each boundary pose.
   *
   * @param[in]  number_of_neighbors  The number of neighbors
   */
  void setNumberOfNeighbors(int number_of_neighbors) { number_of_neighbors_ = number_of_neighbors; }

  /**
   * @brief      Gets the number of neighbors to find for each boundary pose.
   *
   * @return     The number of neighbors.
   */
  int getNumberOfNeighbors(void) { return number_of_neighbors_; }

  /**
   * @brief      Sets the debug display this requires the user to also set the visual cloud.
   *
   * @param[in]  debug_display  The debug display
   */
  void setDebugDisplay(bool debug_display) { if (debug_display) { debug_display_ = true; } }
  
  /**
   * @brief      Sets the visual cloud if debug display is enabled.
   *
   * @param[in]  colored_cloud_ptr  The colored cloud pointer
   */
  void setVisualCloud(const pcl::PointCloud<pcl::PointXYZRGB>::Ptr colored_cloud_ptr)
  {
    visual_cloud_->clear();
    for (const auto &pt : *colored_cloud_ptr) { visual_cloud_->push_back(pt); }
  }
  
private:
  /**
   * @brief      Determines if it an Eigen Matrix4f contains any NaNs
   *
   * @param[in]  matrix  The matrix
   *
   * @return     True if the matrix contains nans, False otherwise.
   */
  static bool containsNaNs(Eigen::Matrix4f matrix);

  /**
   * @brief      Iterates through a vector of poses and remove all NaN's
   *
   * @param[in]  original_boundary_poses  The original boundary poses
   * @param      boundary_poses_no_nan    The boundary poses without NaN's
   */
  static void removeNaNFromPoseTrajectory(const EigenPoseMatrix &original_boundary_poses,
                                          EigenPoseMatrix &boundary_poses_no_nan);

  /**
   * @brief      Iterates through a vector of boundary poses and creates a point cloud at each pose 
   *             of the nearest N points.
   *
   * @param[in]  input_cloud             The input cloud
   * @param[in]  boundary_poses          The boundary poses
   * @param[in]  number_of_neighbors     The number of neighbors
   * @param      boundary_pose_neighbor  Vector of point clouds containing neighbor points for every pose
   */
  static void nearestNNeighborSearch(const pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud,
                                     const EigenPoseMatrix &boundary_poses,
                                     const int &number_of_neighbors,
                                     PointCloudVector &boundary_pose_neighbor);

  /**
   * @brief      Iterates through a vector of boundary poses and creates a point cloud of each pose 
   *             by finding the nearest points within a search radius.
   *             Note: Not currently being used, kept for future reference.
   *
   * @param[in]  input_cloud              The input cloud
   * @param[in]  boundary_poses           The boundary poses
   * @param[in]  search_radius            The search radius
   * @param      boundary_pose_neighbors  Vector of point clouds containing neighbor points for every pose
   */
  static void nearestNeighborRadiusSearch(const pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud,
                                          const EigenPoseMatrix &boundary_poses,
                                          const float &search_radius,
                                          PointCloudVector &boundary_pose_neighbors);

  /**
   * @brief      Calculates the standard deviation given a vector of deviations.
   *
   * @param[in]  deviations  Vector of the deviations of a point from a plane.
   *
   * @return     The allowed deviation.
   */
  static float calculateAllowedDeviation(const std::vector<float> &deviations);

  /**
   * @brief      Calculates the result of a plane if the plane is: a*x + b*y + c*z = d
   *
   * @param[in]  x     x - point
   * @param[in]  y     y - point
   * @param[in]  z     z - point
   * @param[in]  a     a - coefficient
   * @param[in]  b     b - coefficient
   * @param[in]  c     c - coefficient
   * @param[in]  d     d - coefficient
   *
   * @return     a*x + b*y + c*z - d = return
   */
  static float calculatePlane(const float &x, const float &y, const float &z, 
                              const float &a, const float &b, const float &c, const float &d);

  /**
   * @brief      Given a vector of boundary poses and the point cloud at each pose,
   *             this function will refine the point clouds by removing the points that
   *             do not lie within the plane of the original pose within some tolerance.
   *             
   * @details    Given a point A = (Ax, Ay, Az), a normal vector at A, N = <Nx, Ny, Nz>,
   *             and a random point on the plane P = (Px, Py, Pz).
   *             
   *             To check if a point is in a plane:
   *             We know that the dot product between the vector AP and N should equal
   *             zero since they are orthogonal.
   *             
   *             dot(P-A, N) = 0 -> dot(P,N) - dot(A,N) = 0 -> dot(P,N) = dot(A,N)
   *             Since A and N are known, substitute values for P to check if it is on the plane.
   *             
   *             This function also calculates the deviation from the plane of all
   *             points in the point cloud. Then it calculates the standard deviation
   *             of the deviations to determine an allowed error.
   *             
   *             If the points in the point cloud fall within this allowed error, the 
   *             point is added into the refined cloud.                        
   *
   * @param[in]  boundary_poses                  The original boundary poses
   * @param[in]  boundary_pose_neighbor          The vector of point clouds within N neighbors
   * @param      refined_boundary_pose_neighbor  The vector of point clouds within plane
   */
  static void refineNeighborPoints(const EigenPoseMatrix &boundary_poses,
                                   const PointCloudVector &boundary_pose_neighbor,
                                   PointCloudVector &refined_boundary_pose_neighbor);

  /**
   * @brief      Calculates the normals for a point cloud.
   *
   * @param[in]  input_cloud  The input cloud
   * @param      normals      The point cloud of normals
   */
  static void computeNormals(const pcl::PointCloud<pcl::PointXYZ>::Ptr &input_cloud, 
                             pcl::PointCloud<pcl::Normal> &normals);

  /**
   * @brief      Calculates the indicies of boundary points for each point cloud
   *             given a vector of point clouds.
   *
   * @param[in]  refined_cloud           The vector of refined point clouds at each pose
   * @param      boundary_search_radius  The search radius for PCL's boundary estimation
   * @param[in]  refined_boundary        The vector of boundary indices for each pose
   */
  static void computeBoundaryForRefinedCloud(const PointCloudVector &refined_cloud,
                                             const float boundary_search_radius,
                                             PointCloudBoundaryVector &refined_boundary);

  /**
   * @brief      { function_description }
   *
   * @param[in]  refined_points_cloud  The refined points cloud
   * @param[in]  boundary_cloud        The boundary cloud
   * @param      boundary_points       The boundary points
   */
  static void extractBoundaryPointsFromPointCloud(const PointCloudVector &refined_points_cloud,
                                           const PointCloudBoundaryVector &boundary_cloud,
                                           PointCloudVector &boundary_points);

  static bool
  checkIfPointsAreEqual(const pcl::PointXYZ &point_a, const pcl::PointXYZ &point_b);
  // {
  //   if (point_a.x == point_b.x && point_a.y == point_b.y && point_a.z == point_b.z)
  //   {
  //     return true;
  //   }
  //   else { return false; }
  // }

  // This only works for things that go in a circle.
  static pcl::PointCloud<pcl::PointXYZ>
  defineOrderForPointCloud(const pcl::PointCloud<pcl::PointXYZ> &point_cloud);
  // {
  //   pcl::PointCloud<pcl::PointXYZ> unordered_point_cloud;
  //   unordered_point_cloud.reserve(point_cloud.size());
  //   pcl::PointCloud<pcl::PointXYZ> ordered_point_cloud;
  //   ordered_point_cloud.reserve(point_cloud.size());

  //   unordered_point_cloud = point_cloud;
  //   pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;



  //   #if 1
  //   // Safe Method
  //   int K = 2;
  //   std::size_t i = 0;
  //   pcl::PointXYZ searchpoint = point_cloud.points[i];
  //   std::vector<int> pointIdxNKNSearch(K);
  //   std::vector<float> pointNKNSquaredDistance(K);
  //   ordered_point_cloud.push_back(searchpoint);



  //   for (std::size_t j = 0; j < point_cloud.points.size() - 1; j++)
  //   {
  //     kdtree.setInputCloud(unordered_point_cloud.makeShared());
  //     if (kdtree.nearestKSearch(searchpoint, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
  //     {
  //       searchpoint = unordered_point_cloud.points[pointIdxNKNSearch[1]];
  //       ordered_point_cloud.push_back(searchpoint);
  //       unordered_point_cloud.points.erase(unordered_point_cloud.begin() + pointIdxNKNSearch[0]);
  //     }
  //   }
  //   #endif



  //   #if 0
  //   boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer ("temp"));
  //   viewer->setBackgroundColor (0, 0, 0);
  //   pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> unordered(point_cloud.makeShared(), 255, 0, 0);
  //   pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> ordered(ordered_point_cloud.makeShared(), 0, 255, 0);
  //   viewer->addPointCloud<pcl::PointXYZ> (point_cloud.makeShared(), unordered, "unordered cloud");
  //   viewer->addPointCloud<pcl::PointXYZ> (ordered_point_cloud.makeShared(), ordered, "ordered cloud");
  //   viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "unordered cloud");
  //   viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "ordered cloud");

  //   for (std::size_t i = 0; i < ordered_point_cloud.points.size(); i++)
  //   {
  //     std::vector<float> color = getRGB(mapIntensity(i, 0, point_cloud.points.size(), 0, 100));
  //     std::string shape_name = "shape_" + std::to_string(i);
  //     viewer->addSphere(ordered_point_cloud.points[i], 1.0, color[0], color[1], color[2], shape_name);
  //   }
  //   // viewer->removeShape("asdf");

  //   while (!viewer->wasStopped ())
  //   {
  //     viewer->spinOnce (100);
  //     boost::this_thread::sleep (boost::posix_time::microseconds (100000));
  //   }
  //   #endif

  //   return ordered_point_cloud;
  // }

  static std::vector<float>
  getRGB(float intensity);
  // {
  //   std::vector<float> rgb_vals;
  //   rgb_vals.reserve(3);
  //   rgb_vals.push_back(((255*intensity)/100)/255);
  //   rgb_vals.push_back(((255*(100-intensity))/100)/255);
  //   rgb_vals.push_back(0);
  //   return rgb_vals;
  // }

  static float 
  mapIntensity(float x, float in_min, float in_max, float out_min, float out_max);
  // {
  //   return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
  // }


public:

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
        //for (std::size_t j = 0; j < pointIdxNKNSearch.size(); j++)
        {
          temp_point.x = extracted_boundary_points[i].points[pointIdxNKNSearch[K-1]].x;
          temp_point.y = extracted_boundary_points[i].points[pointIdxNKNSearch[K-1]].y;
          temp_point.z = extracted_boundary_points[i].points[pointIdxNKNSearch[K-1]].z;
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

    if (debug_display_data->rendered_additional_shapes_ == true)
    {
      for (std::size_t i = 0; i < debug_display_data->rendered_shape_count_; i++)
      {
        std::string additional_name = "additional_pose_" + std::to_string(i);
        debug_display_data->viewer_->removeShape(additional_name);
      }
      debug_display_data->rendered_additional_shapes_ = false;
    }
          
    pcl::PointXYZ pose_point;
    pose_point.x = debug_display_data->boundary_poses_[debug_display_data->current_pose_index_](0, 3);
    pose_point.y = debug_display_data->boundary_poses_[debug_display_data->current_pose_index_](1, 3);
    pose_point.z = debug_display_data->boundary_poses_[debug_display_data->current_pose_index_](2, 3);
    // Pose Point
    debug_display_data->viewer_->addSphere(pose_point, 0.001*2.5, 1.0, 0.0, 0.0, "pose point");
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
    debug_display_data->viewer_->addSphere(debug_display_data->new_pose_points_[debug_display_data->current_pose_index_], 0.001*2.5, 0.0, 1.0, 0.0, "new point");

    for (std::map<int, PointVector>::const_iterator it = debug_display_data->additional_poses_.begin();
         it != debug_display_data->additional_poses_.end(); it++)
    {  
      if (it->first == debug_display_data->current_pose_index_)// || (it->first + 1) == debug_display_data->current_pose_index_)
      {
        for (std::size_t i = 0; i < it->second.size(); i++)
        {
          std::string additional_name = "additional_pose_" + std::to_string(i);
          // Might be able to remove this.
          if (i == it->second.size()-1)
            debug_display_data->viewer_->addSphere(it->second[i], 0.001*2.5, 0.0, 1.0, 0.0, additional_name);
          else
            debug_display_data->viewer_->addSphere(it->second[i], 0.001*2.5, 0.0, 1.0, 0.0, additional_name);
        }
        debug_display_data->rendered_shape_count_ = it->second.size();
        debug_display_data->rendered_additional_shapes_ = true;
      }
    }
  }

  void
  debugDisplay(const EigenPoseMatrix &boundary_poses,
               const PointCloudVector &boundary_pose_neighbor,
               const PointCloudVector &refined_boundary_pose_neighbor,
               const PointCloudVector &neighbor_boundary_points,
               const PointVector &new_pose_points,
               const EigenPoseMatrix &refined_poses,
               const std::map<int, PointVector> &additional_poses)
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
                                        neighbor_boundary_points, new_pose_points, additional_poses);

    viewer->registerKeyboardCallback(keyboardEventOccurred, static_cast<void *>(&debug_display_data));

    while (!viewer->wasStopped())
    {
      viewer->spinOnce (100);
      boost::this_thread::sleep(boost::posix_time::microseconds (100000));
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
    return round(1.5*(number_of_points - 2)); // -2 because of the original two points.
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

    if (kdtree.nearestKSearch(neighbor_new_pose_points[index], K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
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
      if ((closest_pose_index - pose_index) < ((boundary_points[index].width / 2)))
      {
        std::cout << "The closest way to get to closest_pose_index is to add" << std::endl;
        if (closest_pose_index > num_poses_required)
        {
          std::cout << "Skipping Points" << std::endl;
          int skip_point = (int)floor(closest_pose_index / num_poses_required);
          for (std::size_t i = pose_index; i < closest_pose_index; i += skip_point)
          {
            additional_points.push_back(boundary_points[index].points[i]);
          }
        }
        else
        {
          std::cout << "Not Skipping Points" << std::endl;
          for (std::size_t i = pose_index; i < closest_pose_index; i++)
          {
            additional_points.push_back(boundary_points[index].points[i]);
          }
        }
      }

      // This means the shortest way to the other pose is to move backwards in the boundary.
      else
      {
        std::cout << "The closest way to get to the closest_pose_index is to subtract" << std::endl;
        //if ((closest_pose_index - ((boundary_points[index].width / 2)) > num_poses_required)
        if ((boundary_points[index].width - 1 - closest_pose_index) > num_poses_required)
        {
          std::cout << "Skipping Points" << std::endl;
          //int skip_point = (int)floor((closest_pose_index - ((boundary_points[index].width - 1) / 2)) / num_poses_required);
          int skip_point = (int)floor((boundary_points[index].width - 1 - closest_pose_index) / num_poses_required);

          /*
          for (std::size_t i = 0; i < (closest_pose_index - ((boundary_points[index].width - 1) / 2)); i+= skip_point)
          {
            if (i == 0) { additional_points.push_back(boundary_points[index].points[i]); }
            else
            {
              int point_index = boundary_points[index].width - i - skip_point;
              additional_points.push_back(boundary_points[index].points[point_index]);
            }
          }
          */
          //for (std::size_t i = (boundary_points[index].width); i >= closest_pose_index; i -= skip_point)
          additional_points.push_back(boundary_points[index].points[pose_index]);
          for (std::size_t i = (boundary_points[index].width - 2); i > closest_pose_index; i -= skip_point)
          {
            additional_points.push_back(boundary_points[index].points[i]);
          }
        }
        else
        {
          std::cout << "Not Skipping Points" << std::endl;
          /*
          for (std::size_t i = 0; i < (closest_pose_index - ((boundary_points[index].width - 1) / 2)); i++)
          {
            if (i == 0) { additional_points.push_back(boundary_points[index].points[i]); }
            else
            {
              int point_index = boundary_points[index].width - i;
              additional_points.push_back(boundary_points[index].points[point_index]);
            }
          }
          */
          for (std::size_t i = closest_pose_index; i < (boundary_points[index].width - 1); i++)
          {
            additional_points.push_back(boundary_points[index].points[i]);
          }
          additional_points.push_back(boundary_points[index].points[pose_index]);
        }
      }
    }
    return additional_points;
  }

  static void
  calculateAdditionalPosesRequiredToFillGaps(const PointCloudVector &boundary_points, 
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

  static EigenPoseMatrix
  convertPointsToEigenMatrix(const EigenPoseMatrix &boundary_poses,
                             const PointVector &points,
                             const int &index)
  {
    EigenPoseMatrix additional_poses;
    additional_poses.reserve(points.size());

    Eigen::Matrix4f temp_pose;
    for (std::size_t i = 0; i < points.size(); i++)
    {
      temp_pose = boundary_poses[index];
      temp_pose(0, 3) = points[i].x;
      temp_pose(1, 3) = points[i].y;
      temp_pose(2, 3) = points[i].z;
      additional_poses.push_back(temp_pose);
    }
    return additional_poses;
  }

  static void
  addAdditionalPosesToRefinedPoses(const EigenPoseMatrix &boundary_poses,
                                   const std::map<int, PointVector> &additional_poses,
                                   EigenPoseMatrix &refined_poses)
  {
    int total_points_to_add = 0;

    std::vector<int> additional_pose_indices;
    additional_pose_indices.reserve(additional_poses.size());

    std::vector<PointVector> additional_pose_points;
    additional_pose_points.reserve(additional_poses.size());
    
    for (std::map<int, PointVector>::const_iterator it = additional_poses.begin(); it != additional_poses.end(); it++)
    {
      additional_pose_indices.push_back(it->first);
      additional_pose_points.push_back(it->second);
      total_points_to_add += it->second.size();
    }

    std::vector<int> new_indices;
    new_indices.reserve(additional_poses.size());

    int count;
    for (std::size_t i = 0; i < additional_pose_indices.size(); i++)
    {
      if (i == 0) 
      { 
        new_indices.push_back(additional_pose_indices[i]);
        count = additional_pose_points[i].size();
      }
      else
      {
        new_indices.push_back(additional_pose_indices[i] + count);
        count += additional_pose_points[i].size();
      }
    }

    for (std::size_t i = 0; i < new_indices.size(); i++)
    {
      EigenPoseMatrix temp_additional_pose_matrix = convertPointsToEigenMatrix(boundary_poses, 
                                                                               additional_pose_points[i], 
                                                                               additional_pose_indices[i]);
      EigenPoseMatrix::iterator it;
      it = refined_poses.begin();
      it = refined_poses.insert(it + new_indices[i] + 1, temp_additional_pose_matrix.begin(), temp_additional_pose_matrix.end());
    }

    // Debug Check
    std::cout << "Size of Refined Pose Matrix: " << refined_poses.size() << std::endl;
    for (std::size_t i = 0; i < new_indices.size(); i++)
    {
      std::cout << "Old Index: " << additional_pose_indices[i] <<
      " Points Added: " << additional_pose_points[i].size() << " New Index: " << new_indices[i] << std::endl;
    }
  }



  static float 
  maxValueOfVector(std::vector<float> &vec)
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
  int number_of_neighbors_;
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr visual_cloud_;

  std::size_t current_pose_index_;
  std::size_t num_poses_;
};

} // namespace godel_scan_tools
#endif // EDGE_REFINEMENT_H