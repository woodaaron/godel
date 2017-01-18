#include <godel_scan_tools/edge_refinement.h>

namespace godel_scan_tools
{
  EigenPoseMatrix
  EdgeRefinement::convertPointsToEigenMatrix(const EigenPoseMatrix &boundary_poses,
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
} // namespace godel_scan_tools