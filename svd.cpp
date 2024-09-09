#include <iostream>
#include <Eigen/Dense>

void svdRegistration(const Eigen::MatrixXd& sourcePoints, const Eigen::MatrixXd& targetPoints, Eigen::Matrix3d& R, Eigen::Vector3d& t)
{
    // ��������
    Eigen::Vector3d centroidSource = sourcePoints.colwise().mean();
    Eigen::Vector3d centroidTarget = targetPoints.colwise().mean();

    // ���Ļ�����
    Eigen::MatrixXd centeredSource = sourcePoints.rowwise() - centroidSource.transpose();
    Eigen::MatrixXd centeredTarget = targetPoints.rowwise() - centroidTarget.transpose();

    // ����Э�������
    Eigen::Matrix3d H = centeredSource.transpose() * centeredTarget;

    // ʹ��SVD��������ֵ�ֽ�
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d U = svd.matrixU();
    Eigen::Matrix3d V = svd.matrixV();

    // ������ת����
    R = V * U.transpose();

    // ����ƽ������
    t = centroidTarget - R * centroidSource;
}

int main()
{
      Eigen::MatrixXd sourcePoints(3, 4);  // Դ���ƣ�ÿ�б�ʾһ���������
    sourcePoints << 1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12;

    Eigen::MatrixXd targetPoints(3, 4);  // Ŀ����ƣ�ÿ�б�ʾһ���������
    targetPoints << 2, 4, 6, 8,
        10, 12, 14, 16,
        18, 20, 22, 24;

    Eigen::Matrix3d R;
    Eigen::Vector3d t;

    svdRegistration(sourcePoints, targetPoints, R, t);

    // ������
    std::cout << "Rotation Matrix: \n" << R << std::endl;
    std::cout << "Translation Vector: \n" << t << std::endl;

    return 0;
}