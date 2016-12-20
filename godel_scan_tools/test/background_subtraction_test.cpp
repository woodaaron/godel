#include "gtest/gtest.h"
#include "godel_scan_tools/background_subtraction.h"
class BGSTestFixture : public ::testing::Test { 
public: 
   BGSTestFixture( ) { 
     distance_threshold = 5.0;

     // Generate random pointcloud data, all values positive and less than 1024 
     pcl::PointCloud<pcl::PointXYZ>::Ptr cloud( new pcl::PointCloud<pcl::PointXYZ>);
     background_cloud_ = cloud;
     background_cloud_->width = 1000;
     background_cloud_->height = 1;
     background_cloud_->points.resize (background_cloud_->width * background_cloud_->height);
     
     for (size_t i = 0; i < background_cloud_->points.size (); ++i){
       background_cloud_->points[i].x = 1024.0f * rand () / (RAND_MAX + 1.0f);
       background_cloud_->points[i].y = 1024.0f * rand () / (RAND_MAX + 1.0f);
       background_cloud_->points[i].z = 1024.0f * rand () / (RAND_MAX + 1.0f);
     }

     pcl::PointCloud<pcl::PointXYZ>::Ptr cld(new pcl::PointCloud<pcl::PointXYZ>);
     foreground_cloud_ = cld;
     foreground_cloud_->width = 1000;
     foreground_cloud_->height = 1;
     foreground_cloud_->points.resize (foreground_cloud_->width * foreground_cloud_->height);
     for (size_t i = 0; i < foreground_cloud_->points.size (); ++i){
       foreground_cloud_->points[i].x = background_cloud_->points[i].x;
       foreground_cloud_->points[i].y = background_cloud_->points[i].y;
       foreground_cloud_->points[i].z = background_cloud_->points[i].z + distance_threshold*1.1;
     }

     pcl::PointCloud<pcl::PointXYZ>::Ptr hcld(new pcl::PointCloud<pcl::PointXYZ>);
     half_cloud_ = hcld;
     half_cloud_->width = 1000;
     half_cloud_->height = 1;
     half_cloud_->points.resize (half_cloud_->width * half_cloud_->height);
     for (size_t i = 0; i < half_cloud_->points.size (); ++i){
       if(i%2){
	 half_cloud_->points[i].x = background_cloud_->points[i].x;
	 half_cloud_->points[i].y = background_cloud_->points[i].y + distance_threshold*1.1;
	 half_cloud_->points[i].z = background_cloud_->points[i].z;
       }
       else{
	 half_cloud_->points[i].x = background_cloud_->points[i].x;
	 half_cloud_->points[i].y = background_cloud_->points[i].y;
	 half_cloud_->points[i].z = background_cloud_->points[i].z;
       }
     }

     BS_.add_background_cloud(background_cloud_);
     BS_.set_distance_threshold(distance_threshold);

   } 
 
   void SetUp( ) { 
       // code here will execute just before the test ensues 
   }
 
   void TearDown( ) { 
       // code here will be called just after the test completes
       // ok to through exceptions from here if need be
   }
 
   ~BGSTestFixture( )  { 
       // cleanup any pending stuff, but no exceptions allowed
     background_cloud_->clear();
     foreground_cloud_->clear();
     half_cloud_->clear();
   }
  pcl::PointCloud<pcl::PointXYZ>::Ptr background_cloud_;
  pcl::PointCloud<pcl::PointXYZ>::Ptr foreground_cloud_;
  pcl::PointCloud<pcl::PointXYZ>::Ptr half_cloud_;
  background_subtraction BS_;
  double distance_threshold;
   // put in any custom data members that you need 
};
// test if background cloud is tested against itself, should be no foreground
TEST_F(BGSTestFixture, no_foreground){
  pcl::PointCloud<pcl::PointXYZ> foreground =   BS_.remove_background(background_cloud_);
  EXPECT_EQ(foreground.points.size(), 0);
}
// test if foreground cloud is all in the foreground
TEST_F(BGSTestFixture, all_foreground){
  pcl::PointCloud<pcl::PointXYZ> foreground =   BS_.remove_background(foreground_cloud_);
  EXPECT_EQ(foreground.points.size(), foreground_cloud_->points.size());
}
// test if background cloud is tested against half foreground and half backgound, it should be half sized
TEST_F(BGSTestFixture, half_foreground){
  pcl::PointCloud<pcl::PointXYZ>foreground =   BS_.remove_background(half_cloud_);
  EXPECT_EQ(foreground.points.size(), foreground_cloud_->points.size()/2);
}

int main(int argc, char **argv){
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
