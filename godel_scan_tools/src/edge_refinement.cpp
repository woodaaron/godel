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
    
    computeBoundaryForRefinedCloud(refined_boundary_pose_neighbor, boundary_search_radius_, neighbor_boundary);

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
                                         const int &number_of_neighbors,
                                         PointCloudVector &boundary_pose_neighbor)             
  {
    boundary_pose_neighbor.resize(boundary_poses.size());

    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(input_cloud);

    #pragma omp parallel for
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

      boundary_pose_neighbor[i] = temp_cloud;
    }
  }

  void 
  EdgeRefinement::nearestNeighborRadiusSearch(const pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud,
                                              const EigenPoseMatrix &boundary_poses,
                                              const float &search_radius,
                                              PointCloudVector &boundary_pose_neighbors)
  {
    boundary_pose_neighbors.resize(boundary_poses.size());

    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(input_cloud);

    #pragma omp parallel for
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

      boundary_pose_neighbors[i] = temp_cloud;
    }
  }

  float 
  EdgeRefinement::calculateAllowedDeviation(const std::vector<float> &deviations)
  {
    float sum = std::accumulate(deviations.begin(), deviations.end(), 0.0);
    float mean = sum / deviations.size();

    std::vector<float> diff(deviations.size());
    std::transform(deviations.begin(), deviations.end(), diff.begin(), [mean](double x) { return x - mean; });
    float sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    float allowed_deviation = std::sqrt(sq_sum / deviations.size());

    return allowed_deviation;
  } 

  float
  EdgeRefinement::calculatePlane(const float &x, const float &y, const float &z, 
                 const float &a, const float &b, const float &c, 
                 const float &d)
  {
    return ((a*x + b*y + c*z) - d);
  }

  void 
  EdgeRefinement::refineNeighborPoints(const EigenPoseMatrix &boundary_poses,
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

  void
  EdgeRefinement::computeNormals(const pcl::PointCloud<pcl::PointXYZ>::Ptr &input_cloud, 
                                 pcl::PointCloud<pcl::Normal> &normals)
  {
    pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> normal_estimation;
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());

    normal_estimation.setInputCloud(input_cloud);
    normal_estimation.setSearchMethod(tree);
    normal_estimation.setKSearch(5);
    normal_estimation.compute(normals);
  }

  void
  EdgeRefinement::computeBoundaryForRefinedCloud(const PointCloudVector &refined_cloud,
                                                 const float boundary_search_radius,
                                                 PointCloudBoundaryVector &refined_boundary)
  {
    refined_boundary.resize(refined_cloud.size());
    
    #pragma omp parallel for
    for (std::size_t i = 0; i < refined_cloud.size(); i++)
    {
      pcl::PointCloud<pcl::Boundary> boundaries;
      pcl::PointCloud<pcl::Normal> normals;

      computeNormals(refined_cloud[i].makeShared(), normals);

      pcl::BoundaryEstimation<pcl::PointXYZ, pcl::Normal, pcl::Boundary> boundary_estimation;
      pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());

      boundary_estimation.setInputCloud(refined_cloud[i].makeShared());
      boundary_estimation.setInputNormals(normals.makeShared());
      boundary_estimation.setRadiusSearch(boundary_search_radius);
      boundary_estimation.setSearchMethod(tree);
      boundary_estimation.setAngleThreshold(90.0 * 3.14 / 180.0);
      boundary_estimation.compute(boundaries);

      refined_boundary[i] = boundaries;
    }
  }

  void
  EdgeRefinement::extractBoundaryPointsFromPointCloud(const PointCloudVector &refined_points_cloud,
                                      const PointCloudBoundaryVector &boundary_cloud,
                                      PointCloudVector &boundary_points)
  {
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
      boundary_points.push_back(defineOrderForPointCloud(temp_cloud));
    }     
  }

  pcl::PointCloud<pcl::PointXYZ>
  EdgeRefinement::defineOrderForPointCloud(const pcl::PointCloud<pcl::PointXYZ> &point_cloud)
  {
    pcl::PointCloud<pcl::PointXYZ> unordered_point_cloud;
    unordered_point_cloud.reserve(point_cloud.size());
    pcl::PointCloud<pcl::PointXYZ> ordered_point_cloud;
    ordered_point_cloud.reserve(point_cloud.size());

    unordered_point_cloud = point_cloud;
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;

    int K = 2;
    std::size_t i = 0;
    pcl::PointXYZ searchpoint = point_cloud.points[i];
    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);
    ordered_point_cloud.push_back(searchpoint);



    for (std::size_t j = 0; j < point_cloud.points.size() - 1; j++)
    {
      kdtree.setInputCloud(unordered_point_cloud.makeShared());
      if (kdtree.nearestKSearch(searchpoint, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
      {
        searchpoint = unordered_point_cloud.points[pointIdxNKNSearch[1]];
        ordered_point_cloud.push_back(searchpoint);
        unordered_point_cloud.points.erase(unordered_point_cloud.begin() + pointIdxNKNSearch[0]);
      }
    }

    #if 0
      boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer ("temp"));
      viewer->setBackgroundColor (0, 0, 0);
      pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> unordered(point_cloud.makeShared(), 255, 0, 0);
      pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> ordered(ordered_point_cloud.makeShared(), 0, 255, 0);
      viewer->addPointCloud<pcl::PointXYZ> (point_cloud.makeShared(), unordered, "unordered cloud");
      viewer->addPointCloud<pcl::PointXYZ> (ordered_point_cloud.makeShared(), ordered, "ordered cloud");
      viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "unordered cloud");
      viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "ordered cloud");

      for (std::size_t i = 0; i < ordered_point_cloud.points.size(); i++)
      {
        std::vector<float> color = getRGB(mapIntensity(i, 0, point_cloud.points.size(), 0, 100));
        std::string shape_name = "shape_" + std::to_string(i);
        viewer->addSphere(ordered_point_cloud.points[i], 1.0, color[0], color[1], color[2], shape_name);
      }

      while (!viewer->wasStopped ())
      {
        viewer->spinOnce (100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
      }
    #endif

    return ordered_point_cloud;
  }
  
  std::vector<float>
  EdgeRefinement::getRGB(float intensity)
  {
    std::vector<float> rgb_vals;
    rgb_vals.reserve(3);
    rgb_vals.push_back(((255*intensity)/100)/255);
    rgb_vals.push_back(((255*(100-intensity))/100)/255);
    rgb_vals.push_back(0);
    return rgb_vals;
  }

  float 
  EdgeRefinement::mapIntensity(float x, float in_min, float in_max, float out_min, float out_max)
  {
    return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
  }



} // namespace godel_scan_tools