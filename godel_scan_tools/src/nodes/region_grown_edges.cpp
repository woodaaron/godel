/*
 * Software License Agreement (BSD License)
 *
 * Point Cloud Library (PCL) - www.pointclouds.org
 * Copyright (c) 2009-2011, Willow Garage, Inc.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following
 *   disclaimer in the documentation and/or other materials provided
 *   with the distribution.
 * * Neither the name of Willow Garage, Inc. nor the names of its
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * $Id:$
 *
 */

// STL
#include <iostream>

// PCL
#include <pcl/filters/filter.h>
#include <pcl/features/normal_3d.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/segmentation/region_growing.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/common/time.h>
#include <pcl/console/parse.h>
#include <pcl/visualization/pcl_visualizer.h>
#include "godel_scan_tools/surface_segmentation.h"
#include "godel_scan_tools/background_subtraction.h"
#include "godel_scan_tools/edge_refinement.h"

int
main (int argc, char** av)
{
  if (argc < 3)
  {
    pcl::console::print_info ("Syntax is: %s <source-pcd-background_file> <source-pcd-with_part_file> [-dump] [-select_segment] [-meters]\n\n", av[0]);
    pcl::console::print_info ("If -dump is provided write the extracted clusters to cluster.dat\n\n");
    return (1);
  }

  pcl::PointCloud<pcl::PointXYZ>::Ptr bg_cloud_ptr (new pcl::PointCloud<pcl::PointXYZ>());
  pcl::PointCloud<pcl::PointXYZ>::Ptr bg_cloud_nonans(new pcl::PointCloud<pcl::PointXYZ>());
  pcl::PointCloud<pcl::PointXYZ>::Ptr bg_cloud_nozeros(new pcl::PointCloud<pcl::PointXYZ>());
  pcl::PointCloud<pcl::PointXYZ>::Ptr hd_cloud(new pcl::PointCloud<pcl::PointXYZ>());
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_ptr (new pcl::PointCloud<pcl::PointXYZ>());
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_no_nans (new pcl::PointCloud<pcl::PointXYZ>());
  pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>());
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_segmented (new pcl::PointCloud<pcl::PointXYZRGB>());
  pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud_with_normals (new pcl::PointCloud<pcl::PointXYZRGBNormal>());

  pcl::PCDWriter writer;
  if (pcl::io::loadPCDFile(av[1], *bg_cloud_ptr)==-1) // load the background cloud
  {
    return -1;
  }
  if (pcl::io::loadPCDFile(av[2], *cloud_ptr)==-1) // load the cloud containing part and background
  {
    return -1;
  }

  pcl::console::print_highlight ("Loaded background %s of size %lu and part %s of size %lu\n", 
				 av[1], bg_cloud_ptr->points.size (),
				 av[2], cloud_ptr->points.size ());


  
  // Remove the nans from both the part scan and the background scan
  cloud_ptr->is_dense = false;
  bg_cloud_ptr->is_dense = false;
  cloud_no_nans->is_dense = false;
  bg_cloud_nonans->is_dense = false;
  bg_cloud_nozeros->is_dense = false;
  hd_cloud->is_dense = false;
  std::vector<int> indices, bg_indices;

  // remove NANS from both clouds
  pcl::console::print_highlight ("Remove NANs from clouds\n");
  pcl::removeNaNFromPointCloud (*cloud_ptr, *cloud_no_nans, indices);
  pcl::removeNaNFromPointCloud (*bg_cloud_ptr, *bg_cloud_nonans, bg_indices);

  // remove the background from the part scan
  pcl::console::print_highlight ("subtract background\n");
  background_subtraction BS(bg_cloud_nonans, 5.0);

  for(int i=0; i<cloud_no_nans->points.size(); i+=(rand()%100)*cloud_no_nans->points.size()*.000005)
  {
    cloud_no_nans->points[i].z +=8.0;
  }
  //  pcl::PointCloud<pcl::PointXYZ>::Ptr part_cloud_ptr = BS.remove_background(cloud_no_nans);

  // sub-sample cloud randomly to increase processing speed for testing
  pcl::PointCloud<pcl::PointXYZ>::Ptr part_cloud_ptr(new pcl::PointCloud<pcl::PointXYZ>());

  for(int i=0;i<bg_cloud_nonans->points.size(); i+=rand()%20)
  {
    pcl::PointXYZ pt = bg_cloud_nonans->points[i];
    if(pt.x==0.0 && pt.y==0.0 && pt.z==0.0)
    {
      // do nothing
    }
    else
    {
      part_cloud_ptr->push_back(pt);
    }
  }
  
  for(int i=0;i<bg_cloud_nonans->points.size(); i++)
  {
    pcl::PointXYZ pt = bg_cloud_nonans->points[i];
    if(pt.x==0.0 && pt.y==0.0 && pt.z==0.0)
    {
      // do nothing
    }
    else
    {
      bg_cloud_nozeros->push_back(pt);
      hd_cloud->push_back(pt);
    }
  }

  // Segment the part into surface regions using a "region growing" scheme
  pcl::console::print_highlight ("segmenting\n");
  surfaceSegmentation SS(part_cloud_ptr); // removes NANs and computes normals

  if(pcl::console::find_switch(argc, av, "-meters"))
  {
    SS.setSearchRadius(.03);
  }
  else
  {
    SS.setSearchRadius(30.0);
  }

  pcl::PointCloud<pcl::PointXYZRGB>::Ptr colored_cloud_ptr(new pcl::PointCloud<pcl::PointXYZRGB>());
  std::vector <pcl::PointIndices> clusters = SS.computeSegments(colored_cloud_ptr);
  pcl::console::print_highlight ("Segmented into %d clusters, colored cloud has %d points\n", clusters.size(), colored_cloud_ptr->points.size());

  // Write the resulting setmented cloud into a pcd file
  writer.write<pcl::PointXYZRGB> ("segmented_part.pcd", *colored_cloud_ptr, false);


  // Select the desired surface segment for processing, currently either the default(largest), or specified via cmd line [-select_segment %d] 
  int selected_segment=-1;
  pcl::console::parse_argument (argc, av, "-select_segment", selected_segment);

  if(selected_segment<0)
  {
    int max_size = 0;
    for(int i=0;i<clusters.size();i++)
    {
      if(clusters[i].indices.size() > max_size)
      {
      	max_size = clusters[i].indices.size();
      	selected_segment = i;
      }
    }
  }

  pcl::console::print_highlight ("select_segment %d has %d points\n", selected_segment, clusters[selected_segment].indices.size());
 
  // extract the selected segment from the part cloud
  pcl::PointCloud<pcl::PointXYZ>::Ptr segmented_surface_ptr(new pcl::PointCloud<pcl::PointXYZ>());

  for(int i=0; i<clusters[selected_segment].indices.size(); i++)
  {
    int index = clusters[selected_segment].indices[i];
    pcl::PointXYZ pt(part_cloud_ptr->points[index]);
    segmented_surface_ptr->points.push_back(pt);
  }

  SS.setInputCloud(segmented_surface_ptr);

  pcl::console::print_highlight ("computing boundary of segment cloud\n");
  pcl::PointCloud<pcl::Boundary>::Ptr    boundary_ptr = SS.getBoundaryCloud();	
  pcl::PointCloud<pcl::PointXYZ>::Ptr boundary_cloud_ptr(new pcl::PointCloud<pcl::PointXYZ>());

  int k=0;

  pcl::IndicesPtr boundary_idx(new std::vector<int>());

  BOOST_FOREACH(pcl::Boundary pt, boundary_ptr->points)
  {
    if(pt.boundary_point)
    {
      boundary_cloud_ptr->points.push_back(segmented_surface_ptr->points[k]);
      boundary_idx->push_back(k);
    }
    k++;
  }

  pcl::console::print_highlight ("cloud has %d boundary points\n", boundary_cloud_ptr->points.size());

  // sort the boundaries
  std::vector< pcl::IndicesPtr > sorted_boundaries;
  int num_boundaries = SS.sortBoundary(boundary_idx, sorted_boundaries);
  int max=0;
  int max_idx=0;

  for(int i=0;i<sorted_boundaries.size();i++)
  {
    pcl::console::print_highlight ("boundary %d has %d points\n", i, sorted_boundaries[i]->size());
    if(sorted_boundaries[i]->size() > max)
    {
      max = sorted_boundaries[i]->size();
      max_idx = i;
    }
  }

  // show surface with boundary line drawn
  boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
  viewer->setBackgroundColor (0, 0, 0);

  pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> boundary_color(boundary_cloud_ptr, 255, 255, 255);
  pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> surface_color(segmented_surface_ptr, 55, 76, 150);
  pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(colored_cloud_ptr);
  //  viewer->addPointCloud<pcl::PointXYZ> (boundary_cloud_ptr, boundary_color, "boundary cloud");
  viewer->addPointCloud<pcl::PointXYZRGB> (colored_cloud_ptr, rgb, "colored cloud");
  //./  viewer->addPointCloud<pcl::PointXYZ> (hd_cloud, surface_color, "HD cloud");
  //  viewer->addPointCloud<pcl::PointXYZ> (segmented_surface_ptr, surface_color, "segmented_surface");
  //  viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "segmented_surface");
  viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "colored cloud");
  //./  viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "HD cloud");
  //  viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "boundary cloud");
  //  viewer->addCoordinateSystem (1.0);
  viewer->initCameraParameters ();

  float red[6], green[6], blue[6];
  int color = 0;
  red[0] = 1.0; green[0] = 0.0;  blue[0] = 0.0; 
  red[1] = 0.0; green[1] =1.0;   blue[1] = 0.0;
  red[2] = 0.0; green[2] =0.0;   blue[2] = 1.0; 
  red[3] = 0.0; green[3] =1.0;   blue[3] = 1.0; 
  red[4] = 1.0; green[4] =0.0;   blue[4] = 1.0; 
  red[5] = 1.0; green[5] =1.0;   blue[5] = 0.0; 

  int q = 0;
  for(int i=0;i<  sorted_boundaries.size();i++)
  {
    for(int j=0; j<sorted_boundaries[i]->size()-1; j++)
    {
      char line_number[255];
      sprintf(line_number,"%03d",q++);
      std::string ls = std::string("line_") + std::string(line_number);

      int idx1 = sorted_boundaries[i]->at(j);
      int idx2 = sorted_boundaries[i]->at(j+1);
      viewer->addLine<pcl::PointXYZ> ( segmented_surface_ptr->points[idx1],
				       segmented_surface_ptr->points[idx2],
				      ls.c_str());
      viewer->setShapeRenderingProperties ( pcl::visualization::PCL_VISUALIZER_COLOR, red[color], green[color], blue[color], ls.c_str());
    }// end for each boundary point
    color = (color+1)%6;
  } // end for each boundary 

  // set up smoothing for trajectory
  std::vector<double> filt_coef;
  filt_coef.push_back(1);
  filt_coef.push_back(2);
  filt_coef.push_back(3);
  filt_coef.push_back(4);
  filt_coef.push_back(5);
  filt_coef.push_back(4);
  filt_coef.push_back(3);
  filt_coef.push_back(2);
  filt_coef.push_back(1);
  SS.setSmoothCoef(filt_coef);	// note, automatically normalizes coefficients for unity gain of filter

  //  std::vector<Pose> pose_trajectory;
  std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f> > pose_trajectory;
  SS.getBoundaryTrajectory(sorted_boundaries, 0, pose_trajectory);
  pcl::console::print_highlight ("pose_trajectory has %d poses\n", pose_trajectory.size());
  
  std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f> > refined_pose_trajectory;
  edgeRefinement EF(bg_cloud_nozeros);
  EF.refineBoundary(pose_trajectory, refined_pose_trajectory);

  q = 0;
  //for(int i=0;i<pose_trajectory.size();i++)
  for(int i = 0; i < refined_pose_trajectory.size(); i++)
  {
    if(q++ % 20 == 0)
    {
      char line_number[255];
      sprintf(line_number,"%03d",q++);
      std::string ls = std::string("pose_") + std::string(line_number);
      Eigen::Affine3f pose(refined_pose_trajectory[i].matrix()); // I think this is where it is breaking.
      viewer->addCoordinateSystem (.030, pose, 0);
    }
  }

  color = (color+2)%6;

  // for(int j=0; j<refined_pose_trajectory.size()-1; j++)
  for(int j=0; j<refined_pose_trajectory.size(); j++) // Removing -1 fixes this segfault.
  {
    char line_number[255];
    sprintf(line_number,"%03d",q++);
    std::string ls = std::string("rline_") + std::string(line_number);
    pcl::PointXYZ p1(refined_pose_trajectory[j](0,3),refined_pose_trajectory[j](1,3),refined_pose_trajectory[j](2,3));
    pcl::PointXYZ p2(refined_pose_trajectory[j+1](0,3),refined_pose_trajectory[j+1](1,3),refined_pose_trajectory[j+1](2,3));
    viewer->addLine<pcl::PointXYZ> ( p1,p2,ls.c_str());
    viewer->setShapeRenderingProperties ( pcl::visualization::PCL_VISUALIZER_COLOR, red[color], green[color], blue[color], ls.c_str());
  }// end for each boundary point
  
  
  if (pcl::console::find_switch (argc, av, "-dump"))
  {
    pcl::console::print_highlight ("Writing clusters to clusters.dat\n");
    std::ofstream clusters_file;
    clusters_file.open ("clusters.dat");

    for (std::size_t i = 0; i < clusters.size (); ++i)
    {
      clusters_file << i << "#" << clusters[i].indices.size () << ": ";
      std::vector<int>::const_iterator pit = clusters[i].indices.begin ();
      clusters_file << *pit;

      for (; pit != clusters[i].indices.end (); ++pit)
      {
        clusters_file << " " << *pit;
        clusters_file << std::endl;
      }
    }
    clusters_file.close ();
  }
  
  while (!viewer->wasStopped ())
  {
    viewer->spinOnce (100);
    boost::this_thread::sleep (boost::posix_time::microseconds (100000));
  }
  return (0);
}
