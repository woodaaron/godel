#include <godel_scan_tools/edge_refinement.h>

namespace godel_scan_tools
{
DebugDisplayData::DebugDisplayData(const std::size_t current_pose_index, const std::size_t num_poses, 
                                   pcl::visualization::PCLVisualizer *viewer,
                                   const EigenPoseMatrix boundary_poses, 
                                   const PointCloudVector boundary_pose_neighbor, 
                                   const PointCloudVector refined_boundary_pose_neighbor, 
                                   const PointCloudVector neighbor_boundary_points,
                                   const PointVector new_pose_points,
                                   const std::map<int, PointVector> additional_poses)
{
  current_pose_index_ = current_pose_index;
  num_poses_ = num_poses;
  viewer_ = viewer;

  boundary_poses_ = boundary_poses;
  boundary_pose_neighbor_ = boundary_pose_neighbor;
  refined_boundary_pose_neighbor_ = refined_boundary_pose_neighbor;
  neighbor_boundary_points_ = neighbor_boundary_points;
  new_pose_points_ = new_pose_points;
  additional_poses_ = additional_poses;
  rendered_additional_shapes_ = 0;
  rendered_shape_count_ = 0;
}

EdgeRefinement::EdgeRefinement(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud):
                               tree_(new pcl::search::KdTree<pcl::PointXYZ>()), 
                               point_density_(0), edge_direction_(0)
{
  tree_->setInputCloud(cloud);
  input_cloud_= pcl::PointCloud<pcl::PointXYZ>::Ptr(cloud);
  visual_cloud_ = pcl::PointCloud<pcl::PointXYZRGB>::Ptr(new pcl::PointCloud<pcl::PointXYZRGB>);
  current_pose_index_ = 0;
  debug_display_ = false;
}

  void 
  EdgeRefinement::refineBoundary(const EigenPoseMatrix &original_boundary_poses, 
                                 EigenPoseMatrix &refined_poses)
  {

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

    nearestNNeighborSearch(input_cloud_, boundary_poses, number_of_neighbors_, boundary_pose_neighbor);

    // 2) Narrow down the radius points at each pose to only lie on the x-y plane of the pose with some error.
    PointCloudVector refined_boundary_pose_radius;
    PointCloudVector refined_boundary_pose_neighbor;

    refined_boundary_pose_radius.reserve(boundary_poses.size());
    refined_boundary_pose_neighbor.reserve(boundary_poses.size());

    refineNeighborPoints(boundary_poses, boundary_pose_neighbor, refined_boundary_pose_neighbor);

    // 3) Find all points that are boundaries.
    PointCloudBoundaryVector radius_boundary;
    PointCloudBoundaryVector neighbor_boundary;

    radius_boundary.reserve(boundary_poses.size());
    neighbor_boundary.reserve(boundary_poses.size());
    
    computeBoundaryForRefinedCloud(refined_boundary_pose_neighbor, neighbor_boundary, boundary_search_radius_);

    PointCloudVector radius_boundary_points;
    PointCloudVector neighbor_boundary_points;

    radius_boundary_points.reserve(boundary_poses.size());
    neighbor_boundary_points.reserve(boundary_poses.size());

    extractBoundaryPointsFromPointCloud(refined_boundary_pose_neighbor, neighbor_boundary, neighbor_boundary_points);

    // 4) Find the boundary point that is closest to the original.
    PointVector radius_new_pose_points;
    PointVector neighbor_new_pose_points;

    radius_new_pose_points.reserve(boundary_poses.size());
    neighbor_new_pose_points.reserve(boundary_poses.size());

    calculateClosestPointInBoundaryToPose(boundary_poses, neighbor_boundary_points, neighbor_new_pose_points);

    // 5) Find any outliers in the new pose points (new points that jump too far).
    std::map<int, int> outlier_index; // Pose Number, Number of Points to Add
    calculateOutliersInNewPosePoints(neighbor_new_pose_points, outlier_index);

    // 6) Determines the boundary points that follow the shortest distance between the poses of the two outliers.
    std::map<int, PointVector> additional_poses;
    calculateAdditionalPosesRequiredToFillGaps(neighbor_boundary_points, neighbor_new_pose_points, outlier_index, additional_poses);

    // 7) Move original boundary pose point to new point while keeping same orientation
    movePoseToNewPoint(boundary_poses, neighbor_new_pose_points, refined_poses);

    if (debug_display_)
    {
      debugDisplay(boundary_poses, boundary_pose_neighbor, refined_boundary_pose_neighbor, 
                   neighbor_boundary_points, neighbor_new_pose_points, refined_poses, additional_poses);
    }

    // 8) Determine the indices at whcih to add the additional poses and add in additional poses to the refined poses
    addAdditionalPosesToRefinedPoses(boundary_poses, additional_poses, refined_poses);
  }

  bool
  EdgeRefinement::containsNaNs(Eigen::Matrix4f matrix)
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

  void 
  EdgeRefinement::removeNaNFromPoseTrajectory(const EigenPoseMatrix &original_boundary_poses,
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

  void 
  EdgeRefinement::nearestNNeighborSearch(const pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud,
                                         const EigenPoseMatrix &boundary_poses,
                                         const int number_of_neighbors,
                                         PointCloudVector &boundary_pose_neighbor)
                         
  {
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(input_cloud);

    for (std::size_t i = 0; i < boundary_poses.size(); i++)
    {
      std::vector<int> pointIdxNKNSearch(number_of_neighbors);
      std::vector<float> pointNKNSquaredDistance(number_of_neighbors);

      pcl::PointXYZ searchpoint;
      searchpoint.x = boundary_poses[i](0, 3);
      searchpoint.y = boundary_poses[i](1, 3);
      searchpoint.z = boundary_poses[i](2, 3);

      pcl::PointCloud<pcl::PointXYZ> temp_cloud;
      if (kdtree.nearestKSearch(searchpoint, number_of_neighbors, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
      {
        for (std::size_t j = 0; j < pointIdxNKNSearch.size(); j++)
        {
          temp_cloud.push_back(input_cloud->points[pointIdxNKNSearch[j]]);
        }
      }

      boundary_pose_neighbor.push_back(temp_cloud);
    }
  }

} // namespace godel_scan_tools