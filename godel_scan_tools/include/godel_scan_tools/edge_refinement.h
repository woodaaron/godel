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
#include "godel_scan_tools/surface_segmentation.h"

typedef Eigen::Matrix<float, 1, 3> NormalVector;
typedef Eigen::Matrix<float, 1, 3> PoseOrigin;

/*
    TODO:
    - Calculate boundary point closest to the pose point, keep orientation the same.
    - Typedef all of the long names.
    - Refactor: Iterate through vectors in main function instead of passing vectors around.
*/

class edgeRefinement
{
public:
  /** @brief constructor 
  *   @param cloud input cloud from which you plan to refine a boundary
  */
  edgeRefinement(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud):
                tree_(new pcl::search::KdTree<pcl::PointXYZ>()), point_density_(0), edge_direction_(0)
  {
    tree_->setInputCloud(cloud);
    input_cloud_= pcl::PointCloud<pcl::PointXYZ>::Ptr(cloud);
    //getPointDensity();
  }

  static void
  nearestNeighborRadiusSearch(const pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud,
                              const std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> &boundary_poses,
                              std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> &boundary_pose_neighbors,
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
                         const std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> &boundary_poses,
                         std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> &boundary_pose_neighbor,
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
  refineNeighborPoints(const std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> 
                                &boundary_poses,
                          const std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> 
                                &boundary_pose_neighbor,
                          std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> 
                                &refined_boundary_pose_neighbor)
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
  computeBoundaryForRefinedCloud(const std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> &refined_cloud,
                                 std::vector<pcl::PointCloud<pcl::Boundary>, Eigen::aligned_allocator<pcl::Boundary>> &refined_boundary)
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

  /*
      TODO: Add some sort of ID system???
  */
  static void
  extractBoundaryPointsFromPointCloud(const std::vector<pcl::PointCloud<pcl::PointXYZ>, 
                                            Eigen::aligned_allocator<pcl::PointXYZ>> &refined_points_cloud,
                                      const std::vector<pcl::PointCloud<pcl::Boundary>, 
                                            Eigen::aligned_allocator<pcl::Boundary>> &boundary_cloud,
                                      std::vector<pcl::PointCloud<pcl::PointXYZ>, 
                                            Eigen::aligned_allocator<pcl::PointXYZ>> &boundary_points)
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
  removeNaNFromPoseTrajectory(const std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> &original_boundary_poses,
                              std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> &boundary_poses_no_nan)
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
  calculateClosestPointInBoundaryToPose(const std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> &boundary_poses,
                                        const std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> &extracted_boundary_points,
                                        std::vector<pcl::PointXYZ> &new_pose_points)
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
  comparePoints(const std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> &boundary_poses,
                const std::vector<pcl::PointXYZ> &radius_new_boundary_points,
                const std::vector<pcl::PointXYZ> &neighbor_new_boundary_points)
  {
    assert(boundary_poses.size() == radius_new_boundary_points.size());
    assert(boundary_poses.size() == neighbor_new_boundary_points.size());

    for (size_t i = 0; i < boundary_poses.size(); i++)
    {
      std::cout << "Old Boundary Point: " << "x: " << boundary_poses[i](0, 3) << ", y: " << boundary_poses[i](1, 3) <<
                ", z: " << boundary_poses[i](2, 3) << std::endl;
      std::cout << "Radius New Boundary Point: " << "x: " << radius_new_boundary_points[i].x << ", y: " << radius_new_boundary_points[i].y <<
                ", z: " << radius_new_boundary_points[i].z << std::endl;
      std::cout << "Neighbor New Boundary Point: " << "x: " << neighbor_new_boundary_points[i].x << ", y: " << neighbor_new_boundary_points[i].y <<
                ", z: " << neighbor_new_boundary_points[i].z << std::endl << std::endl;                
    }
  }

  static void
  movePoseToNewPoint(const std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> &boundary_poses,
                     const std::vector<pcl::PointXYZ> &new_boundary_points,
                     std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> &refined_poses)
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
  refineBoundary(const std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> &original_boundary_poses, 
                 std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> &refined_poses)
  {
    refined_poses.clear();

    float search_radius = 30.0;
    int number_of_neighbors = 300;

    // Remove NaN from input boundary poses.
    std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f>> boundary_poses;
    boundary_poses.reserve(original_boundary_poses.size());
    removeNaNFromPoseTrajectory(original_boundary_poses, boundary_poses);

    // 1) Find all points within R1 of each boundary pose.
    std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> boundary_pose_radius;
    std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> boundary_pose_neighbor;
    boundary_pose_radius.reserve(boundary_poses.size());
    boundary_pose_neighbor.reserve(boundary_poses.size());

    nearestNeighborRadiusSearch(input_cloud_, boundary_poses, boundary_pose_radius, search_radius);
    nearestNNeighborSearch(input_cloud_, boundary_poses, boundary_pose_neighbor, number_of_neighbors);

    // 2) Narrow down the radius points at each pose to only lie on the x-y plane of the pose with some error.
    std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> refined_boundary_pose_radius;
    std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> refined_boundary_pose_neighbor;

    refined_boundary_pose_radius.reserve(boundary_poses.size());
    refined_boundary_pose_neighbor.reserve(boundary_poses.size());

    refineNeighborPoints(boundary_poses, boundary_pose_radius, refined_boundary_pose_radius);
    refineNeighborPoints(boundary_poses, boundary_pose_neighbor, refined_boundary_pose_neighbor);

    // 3) Find all points that are boundaries.
    std::vector<pcl::PointCloud<pcl::Boundary>, Eigen::aligned_allocator<pcl::Boundary>> radius_boundary;
    std::vector<pcl::PointCloud<pcl::Boundary>, Eigen::aligned_allocator<pcl::Boundary>> neighbor_boundary;

    radius_boundary.reserve(boundary_poses.size());
    neighbor_boundary.reserve(boundary_poses.size());
    
    computeBoundaryForRefinedCloud(refined_boundary_pose_radius, radius_boundary);
    computeBoundaryForRefinedCloud(refined_boundary_pose_neighbor, neighbor_boundary);

    std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> radius_boundary_points;
    std::vector<pcl::PointCloud<pcl::PointXYZ>, Eigen::aligned_allocator<pcl::PointXYZ>> neighbor_boundary_points;

    radius_boundary_points.reserve(boundary_poses.size());
    neighbor_boundary_points.reserve(boundary_poses.size());

    extractBoundaryPointsFromPointCloud(refined_boundary_pose_radius, radius_boundary, radius_boundary_points);
    extractBoundaryPointsFromPointCloud(refined_boundary_pose_neighbor, neighbor_boundary, neighbor_boundary_points);

    // 4) Find the boundary point that is closest to the original.
    std::vector<pcl::PointXYZ> radius_new_pose_points;
    std::vector<pcl::PointXYZ> neighbor_new_pose_points;

    radius_new_pose_points.reserve(boundary_poses.size());
    neighbor_new_pose_points.reserve(boundary_poses.size());

    calculateClosestPointInBoundaryToPose(boundary_poses, radius_boundary_points, radius_new_pose_points);
    calculateClosestPointInBoundaryToPose(boundary_poses, neighbor_boundary_points, neighbor_new_pose_points);

    // 5) Move original boundary pose point to new point while keeping same orientation
    movePoseToNewPoint(boundary_poses, radius_new_pose_points, refined_poses);

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
    //comparePoints(boundary_poses, radius_new_pose_points, neighbor_new_pose_points);

    //int idx = 2;

    // ????????
    // printf("boundary pose:\n");
    // printf("%7.3f %7.3f %7.3f %7.3f\n", boundary_poses[idx](0,0),boundary_poses[idx](0,1),boundary_poses[idx](0,2), boundary_poses[idx](0,3));
    // printf("%7.3f %7.3f %7.3f %7.3f\n", boundary_poses[idx](1,0),boundary_poses[idx](1,1),boundary_poses[idx](1,2), boundary_poses[idx](1,3));
    // printf("%7.3f %7.3f %7.3f %7.3f\n", boundary_poses[idx](2,0),boundary_poses[idx](2,1),boundary_poses[idx](2,2), boundary_poses[idx](2,3));
    // printf("%7.3f %7.3f %7.3f %7.3f\n", boundary_poses[idx](3,0),boundary_poses[idx](3,1),boundary_poses[idx](3,2), boundary_poses[idx](3,3));
    
    #if 0
    std::cout << "Boundary Pose: " << std::endl;
    std::cout << boundary_poses[idx] << std::endl;

    computeEdgeDirection(boundary_poses[idx]);

    if(edge_direction_ == 0)
    {
      printf("Can't determine edge direction\n");
      return;
    }

    for(int i=0;i<boundary_poses.size();i++)
    {
      Eigen::Matrix4f pose = boundary_poses[i];
      Eigen::Matrix4f rpose;
      refineEdgePose(pose, rpose);
      refined_poses.push_back(rpose);
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

#endif

  #if 0
  /** @brief estimate the number of points per unit volume 
  *   @return number of points per unit vol, units are the same as the inherent distance in the point cloud data
  */
  float 
  getPointDensity()
  {
    if(point_density_>0) return(point_density_);

    int n = input_cloud_->points.size();
    int K=100;

    pcl::PointXYZ pt = input_cloud_->points[n/1];
    std::vector<int> pt_indices(K);
    std::vector<float> pt_sq_distances(K);

    int num_found;
    if((num_found = tree_->nearestKSearch(pt, K, pt_indices, pt_sq_distances))>0)
    {
      double maxd=0;
      int maxi=0;
      for(int i=0;i<K;i++)
      {// note, there should be one point with zero distance
        if(maxd<pt_sq_distances[i]) 
        {
          maxd = pt_sq_distances[i];
          maxi = i;
        }
      }
      printf("maxd = %lf\n",maxd);
      double r = sqrt(maxd);
      double v = 4/3*3.14*r*r*r; /* volume containing K points Kpts/V */
      point_density_ = K/v; // k pts per vol
      radius_ = cbrt(3/(4*3.13*point_density_))/150;
      printf("calculated radius_=%f\n",radius_);
      radius_ = 0.7;
      sradius_ = radius_;
    }
    else
    {
      printf("Could not find %d points near center of input_cloud_\n",K);
      point_density_ = 0.0;
    }

    return(point_density_);
  }

  /** @brief determine if the edge is along positive or negative y direction
  *   @return 1 if along positive y, -1 if along negative y, 0 if not found
  */
  int 
  getEdgeDirection(){return edge_direction_;}

  int 
  computeEdgeDirection(Eigen::Matrix4f near_edge_pose)
  {
    if(getPointDensity()==0)
    {
      printf("without a point density, edge direction cannot be found\n");
      return(0.0);
    }
    else
    {
      printf("point density = %f radius_ =%f\n", point_density_,radius_);
    }

    pcl::PointXYZ pt;
    pt.x = near_edge_pose(0,3);
    pt.y = near_edge_pose(1,3);
    pt.z = near_edge_pose(2,3);

    pcl::PointXYZ yvec;// the y vector should lie on surface in vicinity of pt
    yvec.x = near_edge_pose(0,1);
    yvec.y = near_edge_pose(1,1);
    yvec.z = near_edge_pose(2,1);

    int num_pos=0;
    int num_neg=0;

    for(int i=0;i<10;i++)
    {
      pcl::PointXYZ ptpt;
      pcl::PointXYZ ntpt;
      ptpt.x = pt.x + i*radius_*yvec.x;
      ptpt.y = pt.y + i*radius_*yvec.y;
      ptpt.z = pt.z + i*radius_*yvec.z;
      ntpt.x = pt.x - i*radius_*yvec.x;
      ntpt.y = pt.y - i*radius_*yvec.y;
      ntpt.z = pt.z - i*radius_*yvec.z;
      std::vector<int>pt_indices(10);
      std::vector<float> pt_distsq(10);

      int num_found;
      if((num_found =tree_->radiusSearch(ptpt, sradius_,pt_indices,pt_distsq))>0) 
      {
        num_pos++;
        pt_indices.clear();
      }
      if(tree_->radiusSearch(ntpt, sradius_,pt_indices,pt_distsq)>0)
      {
        num_neg++;
        pt_indices.clear();
      }
    }

    printf("EdgeDirection: num_pos = %d num_neg = %d\n",num_pos,num_neg);

    if(num_pos>num_neg) edge_direction_=-1;
    if(num_neg>num_pos) edge_direction_=1;

    return(edge_direction_);
  }

  void 
  refineEdgePose(Eigen::Matrix4f &near_edge_pose, Eigen::Matrix4f &refined_edge_pose)
  {
    pcl::PointXYZ pt(near_edge_pose(0,3),near_edge_pose(1,3),near_edge_pose(2,3));
    refined_edge_pose = near_edge_pose; // we will only adjust the position

    if(getEdgeDirection()== 0.0)
    {
      computeEdgeDirection(near_edge_pose);

      if(getEdgeDirection() == 0)
      {
        printf("could not ascertain edge direction in vicinity of this pose\n");
        return;
      }
    }

    pcl::PointXYZ yvec(near_edge_pose(0,1),near_edge_pose(1,1),near_edge_pose(0,1));// the y vector lies on surface 
    pcl::PointXYZ tpt; // test point

    int last_index=0; // index to last point found
    int num_out=0;

    tpt.x = pt.x + edge_direction_*yvec.x*num_out*radius_;
    tpt.y = pt.y + edge_direction_*yvec.y*num_out*radius_;
    tpt.z = pt.z + edge_direction_*yvec.z*num_out*radius_;

    std::vector<int> pt_indices;
    std::vector<float> pt_distsq;

    int num_found=0;
    while(((num_found = tree_->radiusSearch(tpt, sradius_,pt_indices,pt_distsq))>0) && num_out<11) // xxxx1xxxxx
    { 
      num_out++;
      last_index = pt_indices[0];
      tpt.x = pt.x + edge_direction_*yvec.x*num_out*radius_;
      tpt.y = pt.y + edge_direction_*yvec.y*num_out*radius_;
      tpt.z = pt.z + edge_direction_*yvec.z*num_out*radius_;
    }

    if(num_out>1)printf("+edge found %d\n ",num_out);

    if(num_out ==11)
    {
      num_out= 1;
      while(((num_found = tree_->radiusSearch(tpt, sradius_,pt_indices,pt_distsq))>0) && num_out<11) // xxxx1xxxxx
      { 
        num_out++;
        last_index = pt_indices[0];
        tpt.x = pt.x - edge_direction_*yvec.x*num_out*radius_;
        tpt.y = pt.y - edge_direction_*yvec.y*num_out*radius_;
        tpt.z = pt.z - edge_direction_*yvec.z*num_out*radius_;
      }
      if(num_out == 11)
      {
        printf("edge not found, no refinement\n");
      }
      else
      {
        printf("-edge found %d\n ",num_out );
      }
    }
    else
    {// update the location to the last point found near the y-axis
      refined_edge_pose(0,3) = input_cloud_->points[last_index].x;
      refined_edge_pose(1,3) = input_cloud_->points[last_index].y;
      refined_edge_pose(2,3) = input_cloud_->points[last_index].z;
    }
  }
  #endif


#if 0
class edgeRefinement
{
public:
  /** @brief constructor 
  *   @param cloud input cloud from which you plan to refine a boundary
  */
  edgeRefinement(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud):
           tree_(new pcl::search::KdTree<pcl::PointXYZ>()), point_density_(0), edge_direction_(0)
  {
    tree_->setInputCloud(cloud);
    input_cloud_= pcl::PointCloud<pcl::PointXYZ>::Ptr(cloud);
    getPointDensity();
  }

  /** @brief estimate the number of points per unit volume 
  *   @return number of points per unit vol, units are the same as the inherent distance in the point cloud data
  */
  float 
  getPointDensity()
  {
    if(point_density_>0) return(point_density_);

    int n = input_cloud_->points.size();
    int K=100;

    pcl::PointXYZ pt = input_cloud_->points[n/1];
    std::vector<int> pt_indices(K);
    std::vector<float> pt_sq_distances(K);

    int num_found;
    if((num_found = tree_->nearestKSearch(pt, K, pt_indices, pt_sq_distances))>0)
    {
      double maxd=0;
      int maxi=0;
      for(int i=0;i<K;i++)
      {// note, there should be one point with zero distance
        if(maxd<pt_sq_distances[i]) 
        {
          maxd = pt_sq_distances[i];
          maxi = i;
        }
      }
      printf("maxd = %lf\n",maxd);
      double r = sqrt(maxd);
      double v = 4/3*3.14*r*r*r; /* volume containing K points Kpts/V */
      point_density_ = K/v; // k pts per vol
      radius_ = cbrt(3/(4*3.13*point_density_))/150;
      printf("calculated radius_=%f\n",radius_);
      radius_ = 0.7;
      sradius_ = radius_;
    }
    else
    {
      printf("Could not find %d points near center of input_cloud_\n",K);
      point_density_ = 0.0;
    }

    return(point_density_);
  }

  /** @brief determine if the edge is along positive or negative y direction
  *   @return 1 if along positive y, -1 if along negative y, 0 if not found
  */
  int 
  getEdgeDirection(){return edge_direction_;}

  int 
  computeEdgeDirection(Eigen::Matrix4f near_edge_pose)
  {
    if(getPointDensity()==0)
    {
      printf("without a point density, edge direction cannot be found\n");
      return(0.0);
    }
    else
    {
      printf("point density = %f radius_ =%f\n", point_density_,radius_);
    }

    pcl::PointXYZ pt;
    pt.x = near_edge_pose(0,3);
    pt.y = near_edge_pose(1,3);
    pt.z = near_edge_pose(2,3);

    pcl::PointXYZ yvec;// the y vector should lie on surface in vicinity of pt
    yvec.x = near_edge_pose(0,1);
    yvec.y = near_edge_pose(1,1);
    yvec.z = near_edge_pose(2,1);

    int num_pos=0;
    int num_neg=0;

    for(int i=0;i<10;i++)
    {
      pcl::PointXYZ ptpt;
      pcl::PointXYZ ntpt;
      ptpt.x = pt.x + i*radius_*yvec.x;
      ptpt.y = pt.y + i*radius_*yvec.y;
      ptpt.z = pt.z + i*radius_*yvec.z;
      ntpt.x = pt.x - i*radius_*yvec.x;
      ntpt.y = pt.y - i*radius_*yvec.y;
      ntpt.z = pt.z - i*radius_*yvec.z;
      std::vector<int>pt_indices(10);
      std::vector<float> pt_distsq(10);

      int num_found;
      if((num_found =tree_->radiusSearch(ptpt, sradius_,pt_indices,pt_distsq))>0) 
      {
        num_pos++;
        pt_indices.clear();
      }
      if(tree_->radiusSearch(ntpt, sradius_,pt_indices,pt_distsq)>0)
      {
        num_neg++;
        pt_indices.clear();
      }
    }

    printf("EdgeDirection: num_pos = %d num_neg = %d\n",num_pos,num_neg);

    if(num_pos>num_neg) edge_direction_=-1;
    if(num_neg>num_pos) edge_direction_=1;

    return(edge_direction_);
  }

  void 
  refineEdgePose(Eigen::Matrix4f &near_edge_pose, Eigen::Matrix4f &refined_edge_pose)
  {
    pcl::PointXYZ pt(near_edge_pose(0,3),near_edge_pose(1,3),near_edge_pose(2,3));
    refined_edge_pose = near_edge_pose; // we will only adjust the position

    if(getEdgeDirection()== 0.0)
    {
      computeEdgeDirection(near_edge_pose);

      if(getEdgeDirection() == 0)
      {
        printf("could not ascertain edge direction in vicinity of this pose\n");
        return;
      }
    }

    pcl::PointXYZ yvec(near_edge_pose(0,1),near_edge_pose(1,1),near_edge_pose(0,1));// the y vector lies on surface 
    pcl::PointXYZ tpt; // test point

    int last_index=0; // index to last point found
    int num_out=0;

    tpt.x = pt.x + edge_direction_*yvec.x*num_out*radius_;
    tpt.y = pt.y + edge_direction_*yvec.y*num_out*radius_;
    tpt.z = pt.z + edge_direction_*yvec.z*num_out*radius_;

    std::vector<int> pt_indices;
    std::vector<float> pt_distsq;

    int num_found=0;
    while(((num_found = tree_->radiusSearch(tpt, sradius_,pt_indices,pt_distsq))>0) && num_out<11) // xxxx1xxxxx
    { 
      num_out++;
      last_index = pt_indices[0];
      tpt.x = pt.x + edge_direction_*yvec.x*num_out*radius_;
      tpt.y = pt.y + edge_direction_*yvec.y*num_out*radius_;
      tpt.z = pt.z + edge_direction_*yvec.z*num_out*radius_;
    }

    if(num_out>1)printf("+edge found %d\n ",num_out);

    if(num_out ==11)
    {
      num_out= 1;
      while(((num_found = tree_->radiusSearch(tpt, sradius_,pt_indices,pt_distsq))>0) && num_out<11) // xxxx1xxxxx
      { 
        num_out++;
        last_index = pt_indices[0];
        tpt.x = pt.x - edge_direction_*yvec.x*num_out*radius_;
        tpt.y = pt.y - edge_direction_*yvec.y*num_out*radius_;
        tpt.z = pt.z - edge_direction_*yvec.z*num_out*radius_;
      }
      if(num_out == 11)
      {
        printf("edge not found, no refinement\n");
      }
      else
      {
        printf("-edge found %d\n ",num_out );
      }
    }
    else
    {// update the location to the last point found near the y-axis
      refined_edge_pose(0,3) = input_cloud_->points[last_index].x;
      refined_edge_pose(1,3) = input_cloud_->points[last_index].y;
      refined_edge_pose(2,3) = input_cloud_->points[last_index].z;
    }
  }

  void refineBoundary(std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f> > &boundary_poses, 
                      std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f> > &refined_poses)
  {
    refined_poses.clear();

    int idx = 2;
    printf("boundary pose:\n");
    printf("%7.3f %7.3f %7.3f %7.3f\n", boundary_poses[idx](0,0),boundary_poses[idx](0,1),boundary_poses[idx](0,2), boundary_poses[idx](0,3));
    printf("%7.3f %7.3f %7.3f %7.3f\n", boundary_poses[idx](1,0),boundary_poses[idx](1,1),boundary_poses[idx](1,2), boundary_poses[idx](1,3));
    printf("%7.3f %7.3f %7.3f %7.3f\n", boundary_poses[idx](2,0),boundary_poses[idx](2,1),boundary_poses[idx](2,2), boundary_poses[idx](2,3));
    printf("%7.3f %7.3f %7.3f %7.3f\n", boundary_poses[idx](3,0),boundary_poses[idx](3,1),boundary_poses[idx](3,2), boundary_poses[idx](3,3));
    
    computeEdgeDirection(boundary_poses[idx]);

    if(edge_direction_ == 0)
    {
      printf("Can't determine edge direction\n");
      return;
    }

    for(int i=0;i<boundary_poses.size();i++)
    {
      Eigen::Matrix4f pose = boundary_poses[i];
      Eigen::Matrix4f rpose;
      refineEdgePose(pose, rpose);
      refined_poses.push_back(rpose);
    }
  }

private:
  pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud_;
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_;
  float point_density_;
  double radius_;
  double sradius_;
  int edge_direction_;
};
#endif