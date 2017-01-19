/*
    TODO:
    - Check ENSENSO scan density and predict the amount of neighbors.
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
   * @brief      Extracts the boundary points from a point cloud given a vector of clouds
   *             containing boundary point indices. Results in a vector of ordered point clouds
   *             of boundary points.
   *
   * @param[in]  refined_points_cloud  The vector of refined point clouds for each pose
   * @param[in]  boundary_cloud        The vector of boundary index clouds for each refined pose
   * @param      boundary_points       The vector of boundary point clouds for each pose
   */
  static void extractBoundaryPointsFromPointCloud(const PointCloudVector &refined_points_cloud,
                                                  const PointCloudBoundaryVector &boundary_cloud,
                                                  PointCloudVector &boundary_points);

  /**
   * @brief      Takes a point cloud of boundary points and orders them. 
   *             
   * @details    This function assumes that these boundary clouds contain points that
   *             form a shape that encloses something, is curved, and is closed.
   *             
   *             The algorithm starts with the default first point of the cloud, and calculates the two
   *             closest points. The first closest point will be itself, and the second closest point will
   *             be the closest point that is not itself. 
   *             
   *             That search point is added to the ordered cloud and removed from the search cloud.
   *             
   *             This is repeated until all the points are removed and we are left with
   *             a ordered point cloud.
   *             
   *             This is used in the future to generate trajectories along the boundaries.
   *             
   *             Note: This function also contains a visualization tool to check if the boundary
   *                   points were actually ordered correctly. This tool is disabled by default.
   *
   * @param[in]  point_cloud  The unordered point cloud
   *
   * @return     The ordered point cloud
   */
  static pcl::PointCloud<pcl::PointXYZ> defineOrderForPointCloud(const pcl::PointCloud<pcl::PointXYZ> &point_cloud);

  /**
   * @brief      Gets a RGB value given an intensity.
   *
   * @param[in]  intensity  The intensity
   *
   * @return     The rgb.
   */
  static std::vector<float> getRGB(float intensity);

  /**
   * @brief      Maps the intensity.
   *
   * @param[in]  x        The value you want to map
   * @param[in]  in_min   In minimum
   * @param[in]  in_max   In maximum
   * @param[in]  out_min  The out minimum
   * @param[in]  out_max  The out maximum
   *
   * @return     The mapped intensity.
   */
  static float mapIntensity(float x, float in_min, float in_max, float out_min, float out_max);

  /**
   * @brief      Calculates the closest point in boundary to pose.
   *
   * @param[in]  boundary_poses             The boundary poses
   * @param[in]  extracted_boundary_points  The extracted boundary points
   * @param      new_pose_points            The new pose points
   */
  static void calculateClosestPointInBoundaryToPose(const EigenPoseMatrix &boundary_poses,
                                                    const PointCloudVector &extracted_boundary_points,
                                                    PointVector &new_pose_points);

  /**
   * @brief      Creates a new pose by moving the position of the original
   *             pose to the closest boundary point while keeping the same orientation.
   *
   * @param[in]  boundary_poses       The original vector of boundary poses
   * @param[in]  new_boundary_points  The vector of closest the closest boundary point to the corresponding pose
   * @param      refined_poses        The refined poses
   */
  static void movePoseToNewPoint(const EigenPoseMatrix &boundary_poses,
                                 const PointVector &new_boundary_points,
                                 EigenPoseMatrix &refined_poses);

  /**
   * @brief      Calculates the index of the refined poses where there is a jump larger than some
   *             prefined threshold.
   *             
   *             NOTE: Debug statements have been left in the implementation since this is still in development.
   *             
   * @details    The outlier index is the first index of the jump. If the index is 100, the large gap is between poses
   *             100 and 101.          
   *
   * @param[in]  neighbor_new_pose_points  The neighbor new pose points
   * @param      outlier_index             The outlier index
   */
  static void calculateOutliersInNewPosePoints(const PointVector &neighbor_new_pose_points,
                                               std::map<int, int> &outlier_index);
  
  /**
   * @brief      Calculates the distance between two points given a point vector and the index of the points.
   *
   * @param[in]  point_vector  The point vector
   * @param[in]  index_1       The index 1
   * @param[in]  index_2       The index 2
   *
   * @return     The distance between two points.
   */
  static float distanceBetweenTwoPoints(const PointVector &point_vector,
                                        const int index_1, const int index_2);

  /**
   * @brief      Calculates the number of points to insert between outliers.
   *
   * @param[in]  distance_between_points  The distance between points
   * @param[in]  standard_deviation       The standard deviation
   *
   * @return     The number of points to insert.
   */
  static int calculateNumberOfPointsToInsert(const float &distance_between_points,
                                             const float &standard_deviation);
 







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