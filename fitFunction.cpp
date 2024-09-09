#include "fitFunction.h"  

//��С���˷����ָ���뾶Բ��ƽ��Բ��
Eigen::VectorXf roundFittingFixedRadius(pcl::PointCloud<pcl::PointXYZ>::ConstPtr cloud, float r)
{
	Eigen::MatrixXf X(cloud->points.size(), 3);
	Eigen::MatrixXf Y(cloud->points.size(), 1);
	for (int i = 0; i < cloud->points.size(); i++)
	{
		X(i, 0) = cloud->points[i].x;
		X(i, 1) = cloud->points[i].y;
		X(i, 2) = 1;

		Y(i, 0) = -pow(cloud->points[i].x, 2) - pow(cloud->points[i].y, 2) - pow(r, 2);
	}
	Eigen::MatrixXf A = (X.transpose() * X).inverse() * X.transpose() * Y;

	//cout << X << endl;

	Eigen::VectorXf CircleCoeff(7);
	CircleCoeff(0) = -A(0, 0) / 2.f;
	CircleCoeff(1) = -A(1, 0) / 2.f;
	CircleCoeff(2) = 0;
	CircleCoeff(3) = 0;
	CircleCoeff(4) = 0;
	CircleCoeff(5) = 1;
	CircleCoeff(6) = r;

	cout << "Բ������Ϊ��(" << CircleCoeff(0) << ", " << CircleCoeff(1) << ")" << endl;

	return CircleCoeff;
}

//��С���˷����Բ��ƽ��Բ��
Eigen::VectorXf roundFitting(pcl::PointCloud<pcl::PointXYZ>::ConstPtr cloud)
{
	Eigen::MatrixXf X(cloud->points.size(), 3);
	Eigen::MatrixXf Y(cloud->points.size(), 1);
	for (int i = 0; i < cloud->points.size(); i++)
	{
		X(i, 0) = cloud->points[i].x;
		X(i, 1) = cloud->points[i].y;
		X(i, 2) = 1;

		Y(i, 0) = -pow(cloud->points[i].x, 2) - pow(cloud->points[i].y, 2);
	}
	Eigen::MatrixXf A = (X.transpose() * X).inverse() * X.transpose() * Y;

	//cout << X << endl;

	Eigen::VectorXf CircleCoeff(7);
	CircleCoeff(0) = -A(0, 0) / 2.f;
	CircleCoeff(1) = -A(1, 0) / 2.f;
	CircleCoeff(2) = 0;
	CircleCoeff(3) = 0;
	CircleCoeff(4) = 0;
	CircleCoeff(5) = 1;
	CircleCoeff(6) = sqrt(-A(2, 0) + pow(CircleCoeff(0), 2) + pow(CircleCoeff(1), 2));

	cout << "Բ������Ϊ��(" << CircleCoeff(0) << ", " << CircleCoeff(1) << ")" << endl;
	cout << "�뾶Ϊ��" << CircleCoeff(6) << endl;

	return CircleCoeff;
}