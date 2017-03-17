#ifndef BACKGROUND_SUBTRACTION_H
#define BACKGROUND_SUBTRACTION_H

#include <pcl_ros/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <boost/foreach.hpp>


/** @class world_background_subtraction
      @brief Maintains record of baseline sensor data to provide method to remove them leaving only new objects in the scene
*/
class background_subtraction{
 public:
  /** @brief default constructor */
  background_subtraction()
    {
      pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
      bg_cloud_ = cloud;
    }

  /** @brief constructor that sets the background cloud, also initializes the KdTree for searching
      @param bg_cloud the set of points defining the background
  */
  background_subtraction(pcl::PointCloud<pcl::PointXYZ>::Ptr bg_cloud, double distance_threshold): distance_threshold_(distance_threshold)
    {
      pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
      bg_cloud_ = cloud;
      set_background_cloud(bg_cloud);
    };

  /** @brief distructor */
  ~background_subtraction()
    {
      bg_cloud_->clear();
    };

  /** @brief sets the background cloud, replaces whatever points exists if any
      @param background_cloud the cloud representing the background
  */
  void set_background_cloud(pcl::PointCloud<pcl::PointXYZ>::Ptr bg_cloud){
    bg_cloud_->clear();
    BOOST_FOREACH(pcl::PointXYZ pt, *bg_cloud){
      bg_cloud_->push_back(pt);
    }
    kd_tree_.setInputCloud(bg_cloud_);
  }

  /** @brief adds new points to the background, and reinitializes the kd_tree for searching
      @param bg_cloud additional background points
  */
  void add_background_cloud(pcl::PointCloud<pcl::PointXYZ>::Ptr bg_cloud){
    BOOST_FOREACH(pcl::PointXYZ pt, *bg_cloud){
      bg_cloud_->push_back(pt);
    }
    kd_tree_.setInputCloud(bg_cloud_);
  }

  /** @brief gets the distance threshold 
      @returns distance_threshold the radious within which a point will be considered part of the background
  */
  double get_distance_threshold(){
    return(distance_threshold_);
  }

  /** @brief sets the distance threshold 
      @param distance_threshold the radious within which a point will be considered part of the background
  */
  void set_distance_threshold(double distance_threshold){
    distance_threshold_ = distance_threshold;
  }

/** @brief uses a kd-tree to extract the foreground from a background point cloud, 
    that being all points further away than the radius defined by distance_threshold_ from all points in bg_cloud_
    @param test_cloud input from which we want the foreground
    @returns the forground points
 */
  pcl::PointCloud<pcl::PointXYZ>::Ptr remove_background(pcl::PointCloud<pcl::PointXYZ>::Ptr test_cloud)
    {
      pcl::PointCloud<pcl::PointXYZ>::Ptr foreground_cloud(new pcl::PointCloud<pcl::PointXYZ>());
      
      double rsquared;
      std::vector<int> close_point_idx;
      std::vector<float> close_point_dist;
      int q=0;
      BOOST_FOREACH(pcl::PointXYZ searchPoint, *test_cloud){
	if ( kd_tree_.radiusSearch (searchPoint, distance_threshold_, close_point_idx, close_point_dist) == 0 )
	  {
	    //pcl::console::print_highlight ("adding a point %d %f %f %f\n", q, searchPoint.x, searchPoint.y, searchPoint.z);
	    foreground_cloud->push_back(searchPoint);
	  }
	q++;
      }
      return foreground_cloud;
    }

 private:
      double distance_threshold_;
      pcl::KdTreeFLANN<pcl::PointXYZ> kd_tree_;
      pcl::PointCloud<pcl::PointXYZ>::Ptr bg_cloud_;
};
#endif
