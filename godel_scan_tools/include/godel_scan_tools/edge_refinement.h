#ifndef EDGE_REFINEMENT_H
#define EDGE_REFINEMENT_H

#include <pcl_ros/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d.h>
#include <pcl/surface/gp3.h>
#include <pcl/features/boundary.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/io/vtk_io.h>
#include <pcl/surface/concave_hull.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <boost/foreach.hpp>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <math.h>
//#include "godel_scan_tools/surface_segmentation.h"

namespace godel_scan_tools
{

typedef std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> PointCloudVector;
typedef std::vector<pcl::PointCloud<pcl::Boundary>, Eigen::aligned_allocator<pcl::Boundary>> PointCloudBoundaryVector;
typedef std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> EigenPoseMatrix;
typedef std::vector<pcl::PointXYZ> PointVector;
typedef Eigen::Matrix<float, 1, 3> NormalVector;
typedef Eigen::Matrix<float, 1, 3> PoseOrigin;

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
  }

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

      float allowed_error = 0.15;

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
                                 PointCloudBoundaryVector &refined_boundary)
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
      boundary_estimation.setRadiusSearch(5.0);
      boundary_estimation.setSearchMethod(tree);
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
  refineBoundary(const EigenPoseMatrix &original_boundary_poses, 
                 EigenPoseMatrix &refined_poses)
  {
    refined_poses.clear();

    float search_radius = 30.0;
    int number_of_neighbors = 500;

    // Remove NaNs from input boundary poses.
    EigenPoseMatrix boundary_poses;
    boundary_poses.reserve(original_boundary_poses.size());
    removeNaNFromPoseTrajectory(original_boundary_poses, boundary_poses);

    // 1) Find all points within R1 of each boundary pose.
    PointCloudVector boundary_pose_radius;
    PointCloudVector boundary_pose_neighbor;
    boundary_pose_radius.reserve(boundary_poses.size());
    boundary_pose_neighbor.reserve(boundary_poses.size());

    //nearestNeighborRadiusSearch(input_cloud_, boundary_poses, boundary_pose_radius, search_radius);
    nearestNNeighborSearch(input_cloud_, boundary_poses, boundary_pose_neighbor, number_of_neighbors);

    // 2) Narrow down the radius points at each pose to only lie on the x-y plane of the pose with some error.
    PointCloudVector refined_boundary_pose_radius;
    PointCloudVector refined_boundary_pose_neighbor;

    refined_boundary_pose_radius.reserve(boundary_poses.size());
    refined_boundary_pose_neighbor.reserve(boundary_poses.size());

    //refineNeighborPoints(boundary_poses, boundary_pose_radius, refined_boundary_pose_radius);
    refineNeighborPoints(boundary_poses, boundary_pose_neighbor, refined_boundary_pose_neighbor);

    // 3) Find all points that are boundaries.
    PointCloudBoundaryVector radius_boundary;
    PointCloudBoundaryVector neighbor_boundary;

    radius_boundary.reserve(boundary_poses.size());
    neighbor_boundary.reserve(boundary_poses.size());
    
    //computeBoundaryForRefinedCloud(refined_boundary_pose_radius, radius_boundary);
    computeBoundaryForRefinedCloud(refined_boundary_pose_neighbor, neighbor_boundary);

    PointCloudVector radius_boundary_points;
    PointCloudVector neighbor_boundary_points;

    radius_boundary_points.reserve(boundary_poses.size());
    neighbor_boundary_points.reserve(boundary_poses.size());

    //extractBoundaryPointsFromPointCloud(refined_boundary_pose_radius, radius_boundary, radius_boundary_points);
    extractBoundaryPointsFromPointCloud(refined_boundary_pose_neighbor, neighbor_boundary, neighbor_boundary_points);

    // 4) Find the boundary point that is closest to the original.
    PointVector radius_new_pose_points;
    PointVector neighbor_new_pose_points;

    radius_new_pose_points.reserve(boundary_poses.size());
    neighbor_new_pose_points.reserve(boundary_poses.size());

    //calculateClosestPointInBoundaryToPose(boundary_poses, radius_boundary_points, radius_new_pose_points);
    calculateClosestPointInBoundaryToPose(boundary_poses, neighbor_boundary_points, neighbor_new_pose_points);

    // 5) Move original boundary pose point to new point while keeping same orientation
    //movePoseToNewPoint(boundary_poses, radius_new_pose_points, refined_poses);
    movePoseToNewPoint(boundary_poses, neighbor_new_pose_points, refined_poses);

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
};
} // namespace godel_scan_tools
#endif // EDGE_REFINEMENT_H