#include <iostream>
#include <Eigen/Dense>

void svdRegistration(const Eigen::MatrixXd& sourcePoints, const Eigen::MatrixXd& targetPoints, Eigen::Matrix3d& R, Eigen::Vector3d& t)
{
    // 计算质心
    Eigen::Vector3d centroidSource = sourcePoints.colwise().mean();
    Eigen::Vector3d centroidTarget = targetPoints.colwise().mean();

    // 中心化点云
    Eigen::MatrixXd centeredSource = sourcePoints.rowwise() - centroidSource.transpose();
    Eigen::MatrixXd centeredTarget = targetPoints.rowwise() - centroidTarget.transpose();

    // 计算协方差矩阵
    Eigen::Matrix3d H = centeredSource.transpose() * centeredTarget;

    // 使用SVD进行奇异值分解
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d U = svd.matrixU();
    Eigen::Matrix3d V = svd.matrixV();

    // 计算旋转矩阵
    R = V * U.transpose();

    // 计算平移向量
    t = centroidTarget - R * centroidSource;
}

int main()
{
      Eigen::MatrixXd sourcePoints(3, 4);  // 源点云，每列表示一个点的坐标
    sourcePoints << 1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12;

    Eigen::MatrixXd targetPoints(3, 4);  // 目标点云，每列表示一个点的坐标
    targetPoints << 2, 4, 6, 8,
        10, 12, 14, 16,
        18, 20, 22, 24;

    Eigen::Matrix3d R;
    Eigen::Vector3d t;

    svdRegistration(sourcePoints, targetPoints, R, t);

    // 输出结果
    std::cout << "Rotation Matrix: \n" << R << std::endl;
    std::cout << "Translation Vector: \n" << t << std::endl;

    return 0;
}