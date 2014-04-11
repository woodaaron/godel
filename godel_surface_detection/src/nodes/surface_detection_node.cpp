/*
	Copyright Feb 12, 2014 Southwest Research Institute

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

		http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
*/

#include <ros/ros.h>
#include <godel_surface_detection/detection/surface_detection.h>
#include <godel_surface_detection/interactive/interactive_surface_server.h>
#include <pcl_ros/transforms.h>
#include <pcl/common/common.h>
#include <pcl/filters/filter.h>
#include <pcl/console/parse.h>

using namespace godel_surface_detection::detection;

const float DEFAULT_ACQUISITION_TIME = 2.0f; //second
const std::string DEFAULT_FRAME_ID = "world_frame";

const std::string POINT_CLOUD_TOPIC="sensor_point_cloud";
const std::string SEGMENTS_CLOUD_TOPIC = "segments_cloud";
const std::string MARKER_ARRAY_TOPIC = "segment_markers";
const std::string NODE_NAME = "surface_detection_node";
const std::string HELP_TEXT="\n" + NODE_NAME + " help:\n" +
		"-h help menu\n" +
		"-a <acquisition time (sec)>\n" +
		"-i <frame id>\n";

// surface detection instance
godel_surface_detection::detection::SurfaceDetection SurfDetect;

// marker server instance
godel_surface_detection::interactive::InteractiveSurfaceServer SurfServer;

// functions
bool load_parameters(std::string node_ns = "")
{
	ros::NodeHandle nh(node_ns.empty() ? "~" : node_ns);
	bool succeeded;
	const std::string ns = params::PARAMETER_NS + "/";
	if(		nh.getParam(ns + params::FRAME_ID,SurfDetect.frame_id_)&&
			nh.getParam(ns + params::STOUTLIER_MEAN,SurfDetect.meanK_) &&
			nh.getParam(ns + params::STOUTLIER_STDEV_THRESHOLD,SurfDetect.stdv_threshold_) &&
			nh.getParam(ns + params::REGION_GROWING_MIN_CLUSTER_SIZE,SurfDetect.rg_min_cluster_size_) &&
			nh.getParam(ns + params::REGION_GROWING_MAX_CLUSTER_SIZE,SurfDetect.rg_max_cluster_size_) &&
			nh.getParam(ns + params::REGION_GROWING_NEIGHBORS,SurfDetect.rg_neightbors_) &&
			nh.getParam(ns + params::REGION_GROWING_SMOOTHNESS_THRESHOLD,SurfDetect.rg_smoothness_threshold_) &&
			nh.getParam(ns + params::REGION_GROWING_CURVATURE_THRESHOLD,SurfDetect.rg_curvature_threshold_) &&
			nh.getParam(ns + params::TRIANGULATION_SEARCH_RADIUS,SurfDetect.tr_search_radius_) &&
			nh.getParam(ns + params::TRIANGULATION_MU ,SurfDetect.tr_mu_) &&
			nh.getParam(ns + params::TRIANGULATION_MAX_NEAREST_NEIGHBORS,SurfDetect.tr_max_nearest_neighbors_) &&
			nh.getParam(ns + params::TRIANGULATION_MAX_SURFACE_ANGLE,SurfDetect.tr_max_surface_angle_) &&
			nh.getParam(ns + params::TRIANGULATION_MIN_ANGLE,SurfDetect.tr_min_angle_) &&
			nh.getParam(ns + params::TRIANGULATION_MAX_ANGLE,SurfDetect.tr_max_angle_) &&
			nh.getParam(ns + params::TRIANGULATION_NORMAL_CONSISTENCY,SurfDetect.tr_normal_consistency_) &&
			nh.getParam(ns + params::VOXEL_LEAF_SIZE,SurfDetect.voxel_leafsize_) &&
			nh.getParam(ns + params::OCCUPANCY_THRESHOLD,SurfDetect.occupancy_threshold_) &&
			nh.getParam(ns + params::MLS_UPSAMPLING_RADIUS,SurfDetect.mls_upsampling_radius_) &&
			nh.getParam(ns + params::MLS_POINT_DENSITY,SurfDetect.mls_point_density_) &&
			nh.getParam(ns + params::MLS_SEARCH_RADIUS,SurfDetect.mls_search_radius_) &&
			nh.getParam(ns + params::USE_TABLETOP_SEGMENTATION,SurfDetect.use_tabletop_seg_) &&
			nh.getParam(ns + params::TABLETOP_SEG_DISTANCE_THRESH,SurfDetect.tabletop_seg_distance_threshold_) &&
			nh.getParam(ns + params::MARKER_ALPHA,SurfDetect.marker_alpha_) &&
			nh.getParam(ns + params::IGNORE_LARGEST_CLUSTER,SurfDetect.ignore_largest_cluster_)
			)
	{
		succeeded = true;
		ROS_INFO_STREAM("surface detection parameters loaded");
	}
	else
	{
		succeeded = false;
		ROS_ERROR_STREAM("surface detection parameter(s) not found");
	}

	return succeeded;
}

void point_cloud_subscriber(const sensor_msgs::PointCloud2ConstPtr msg)
{
	static tf::TransformListener tf_listener;

	// convert to message to point cloud
	Cloud::Ptr new_cloud_ptr(new Cloud());
	pcl::fromROSMsg<pcl::PointXYZ>(*msg,*new_cloud_ptr);

	// removed nans
	std::vector<int> index;
	pcl::removeNaNFromPointCloud(*new_cloud_ptr,*new_cloud_ptr,index);

	// transform to frame id
	tf::StampedTransform source_to_target_tf;
	if(SurfDetect.frame_id_.compare(msg->header.frame_id) !=0)
	{
		ROS_INFO_STREAM("Source cloud with frame id '"<<msg->header.frame_id<<"' will be transformed to frame id: '"
				<<SurfDetect.frame_id_<<"'");
		try
		{
			tf_listener.lookupTransform(SurfDetect.frame_id_,msg->header.frame_id,
					ros::Time::now() - ros::Duration(0.2f),source_to_target_tf);
			pcl_ros::transformPointCloud(*new_cloud_ptr,*new_cloud_ptr,source_to_target_tf);
		}
		catch(tf::LookupException &e)
		{
			ROS_ERROR_STREAM("Transform lookup error, using source frame id '"<< msg->header.frame_id<<"'");
			SurfDetect.frame_id_ = msg->header.frame_id;
		}
		catch(tf::ExtrapolationException &e)
		{
			ROS_ERROR_STREAM("Transform lookup error, using source frame id '"<< msg->header.frame_id<<"'");
			SurfDetect.frame_id_ = msg->header.frame_id;
		}

	}
	else
	{
		ROS_INFO_STREAM("Source cloud is already in frame id '"<<msg->header.frame_id<<", skipping transform");
	}

	// add point cloud
	SurfDetect.add_cloud(*new_cloud_ptr);
}

int main(int argc,char** argv)
{
	ros::init(argc,argv,"surface_detection_node");
	ros::NodeHandle nh;

	// parsing arguments
	float acquisition_time = DEFAULT_ACQUISITION_TIME;
	std::string frame_id = DEFAULT_FRAME_ID;

	if(pcl::console::find_switch(argc,argv,"-h"))
	{
		std::cout<<HELP_TEXT<<std::endl;
		return 0;
	}

	if(pcl::console::find_switch(argc,argv,"-a"))
	{
		pcl::console::parse(argc,argv,"-a",acquisition_time);
	}

	if(pcl::console::find_switch(argc,argv,"-i"))
	{
		pcl::console::parse(argc,argv,"-i",frame_id);
	}

	ROS_INFO_STREAM("Using acquisition time '"<<acquisition_time<<"'");

	// publishers
	ros::Publisher point_cloud_publisher = nh.advertise<sensor_msgs::PointCloud2>(
			SEGMENTS_CLOUD_TOPIC,1);

	ros::Publisher markers_publisher = nh.advertise<visualization_msgs::MarkerArray>(
			MARKER_ARRAY_TOPIC,1);


	// load paramters
	if(!load_parameters() || !SurfDetect.init() || !SurfServer.init())
	{
		ROS_ERROR_STREAM("Initialization failed");
		return 0;
	}

	// start server
	SurfServer.run();

	// acquire data
	ros::Duration total_duration(acquisition_time);
	ros::Time start_time = ros::Time::now();
	ros::Subscriber cloud_subs = nh.subscribe(POINT_CLOUD_TOPIC,1,point_cloud_subscriber);

	ROS_INFO_STREAM("Started point cloud acquisition");
	while(ros::ok() &&
			((start_time + total_duration) > ros::Time::now()))
	{
		ros::spinOnce();
	}
	cloud_subs.shutdown();

	// processing data
	bool succeeded = SurfDetect.acquired_clouds_counter_ >  0;
	if(succeeded)
	{
		if(SurfDetect.find_surfaces())
		{
			ROS_INFO_STREAM(SurfDetect.get_results_summary());

		}
		else
		{
			ROS_ERROR_STREAM("No surface was found");
			succeeded = false;
		}
	}
	else
	{
		ROS_ERROR_STREAM("Data acquisition failed");
		succeeded = false;
	}

	if(succeeded)
	{
		ROS_INFO_STREAM("Publishing segments visuals");
		sensor_msgs::PointCloud2 cloud_msg;
		visualization_msgs::MarkerArray markers_msg = SurfDetect.get_surface_markers();
		SurfDetect.get_region_colored_cloud(cloud_msg);

		// adding markers to server
		for(int i =0;i < markers_msg.markers.size();i++)
		{
			SurfServer.add_surface(markers_msg.markers[i]);
		}

		ros::Duration loop_rate(0.5f);
		while(succeeded && ros::ok() )
		{
			point_cloud_publisher.publish(cloud_msg);
			markers_publisher.publish(markers_msg);
			ros::spinOnce();

			loop_rate.sleep();
		}

		point_cloud_publisher.shutdown();
		markers_publisher.shutdown();
	}

	return 0;

}




