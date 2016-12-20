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
#include <boost/foreach.hpp>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <math.h>

class edgeRefinement{
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
    float getPointDensity()
    {
      if(point_density_>0) return(point_density_);
      int n = input_cloud_->points.size();
      int K=100;
      pcl::PointXYZ pt = input_cloud_->points[n/1];
      std::vector<int> pt_indices(K);
      std::vector<float> pt_sq_distances(K);
      int num_found;
      if((num_found = tree_->nearestKSearch(pt, K, pt_indices, pt_sq_distances))>0){
	double maxd=0;
	int maxi=0;
	for(int i=0;i<K;i++){// note, there should be one point with zero distance
	  if(maxd<pt_sq_distances[i]) {
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
      else{
	printf("Could not find %d points near center of input_cloud_\n",K);
	point_density_ = 0.0;
      }
      return(point_density_);
    }
    /** @brief determine if the edge is along positive or negative y direction
     *   @return 1 if along positive y, -1 if along negative y, 0 if not found
     */
    int getEdgeDirection(){return edge_direction_;}
    int computeEdgeDirection(Eigen::Matrix4f near_edge_pose)
    {
      if(getPointDensity()==0){
	printf("without a point density, edge direction cannot be found\n");
	return(0.0);
      }else{
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
      for(int i=0;i<10;i++){
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
	if((num_found =tree_->radiusSearch(ptpt, sradius_,pt_indices,pt_distsq))>0) {
	  num_pos++;
	  pt_indices.clear();
	}
	if(tree_->radiusSearch(ntpt, sradius_,pt_indices,pt_distsq)>0){
	  num_neg++;
	  pt_indices.clear();
	}
      }
      printf("EdgeDirection: num_pos = %d num_neg = %d\n",num_pos,num_neg);
      if(num_pos>num_neg) edge_direction_=-1;
      if(num_neg>num_pos) edge_direction_=1;
      return(edge_direction_);
    }
      
    void refineEdgePose(Eigen::Matrix4f &near_edge_pose, Eigen::Matrix4f &refined_edge_pose)
    {
      pcl::PointXYZ pt(near_edge_pose(0,3),near_edge_pose(1,3),near_edge_pose(2,3));
      refined_edge_pose = near_edge_pose; // we will only adjust the position
      if(getEdgeDirection()==0.0){
	computeEdgeDirection(near_edge_pose);
	if(getEdgeDirection() ==0){
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
      while(((num_found = tree_->radiusSearch(tpt, sradius_,pt_indices,pt_distsq))>0) && num_out<11){ // xxxx1xxxxx
	num_out++;
	last_index = pt_indices[0];
	tpt.x = pt.x + edge_direction_*yvec.x*num_out*radius_;
	tpt.y = pt.y + edge_direction_*yvec.y*num_out*radius_;
	tpt.z = pt.z + edge_direction_*yvec.z*num_out*radius_;
      }
      if(num_out>1)printf("+edge found %d\n ",num_out );
      if(num_out ==11){
	num_out= 1;
	while(((num_found = tree_->radiusSearch(tpt, sradius_,pt_indices,pt_distsq))>0) && num_out<11){ // xxxx1xxxxx
	  num_out++;
	  last_index = pt_indices[0];
	  tpt.x = pt.x - edge_direction_*yvec.x*num_out*radius_;
	  tpt.y = pt.y - edge_direction_*yvec.y*num_out*radius_;
	  tpt.z = pt.z - edge_direction_*yvec.z*num_out*radius_;
	}
	if(num_out ==11){
	  printf("edge not found, no refinement\n");
	}
	else{
	  printf("-edge found %d\n ",num_out );
	}
      }else{// update the location to the last point found near the y-axis
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
      if(edge_direction_ == 0){
	printf("Can't determine edge direction\n");
	return;
      }
      for(int i=0;i<boundary_poses.size();i++){
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
