#define _CRT_SECURE_NO_WARNINGS 1 

#include "4.0Function.h"
#include "calibFunction.h"
#include "fitFunction.h"  
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <cassert>
#include <string>
#include <cstdlib>
#include <random>
#include <cmath>


constexpr int f() noexcept { return 0; }
using namespace std;
using namespace cv;


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


void addGaussianNoise(Eigen::MatrixXd& points, double stddev) {
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0, stddev);

	for (int i = 0; i < points.cols(); ++i) {
		for (int j = 0; j < points.rows(); ++j) {
			points(j, i) += distribution(generator);
		}
	}
}

int main()
{
	vector<SLSLens> lens;
	getLens(lens);

	//��ͼ
	//string filenameL = "globalMarkerImgL.png";
	//string filenameR = "globalMarkerImgR.png";
	//string filenameLR = "globalMarkerImgLR.png";
	//string filenameRR = "globalMarkerImgRR.png";

	// filenameL = "L1001.BMP";
	//string filenameR = "R1001.BMP";

	string filenameL = "rightL10.BMP";
	string filenameR = "rightR10.BMP";
	string filenameLR = "1L0.BMP";
	string filenameRR = "1R0.BMP";

	cv::Mat imageL = cv::imread(filenameL, 0);
	cv::Mat imageR = cv::imread(filenameR, 0);
	cv::Mat rimageL = cv::imread(filenameLR, 0);
	cv::Mat rimageR = cv::imread(filenameRR, 0);


	// �ؽ���ά��
	vector<glm::dvec3> rPoints;
	vector<glm::dvec3> lPoints;
	vector<glm::dvec2> matchMarkersL, matchMarkersR;
	vector<glm::dvec2> matchMarkersLR, matchMarkersRR;
	vector<double> epipL, epipR;
	vector<double> epipLR, epipRR;
	vector<glm::dvec4>Transformation;

	lPoints = reconstruct3dPoints(lens, imageL, imageR,
		matchMarkersL, matchMarkersR, epipL, epipR);        //ת�����

	rPoints = reconstruct3dPoints(lens, rimageL, rimageR,
		matchMarkersL, matchMarkersRR, epipLR, epipRR);  //��ӦB���������ά��


	//---------------------------------------------------------��ת������---------------------------------------------------------

	Eigen::MatrixXd emarkerPoints(lPoints.size(), 3);
	for (size_t i = 0; i < lPoints.size(); i++)
	{
		emarkerPoints.row(i) << lPoints[i].x, lPoints[i].y, lPoints[i].z;
	}

	Eigen::MatrixXd esourcePoints(rPoints.size(), 3);    
	for (size_t i = 0; i < rPoints.size(); i++)
	{
		esourcePoints.row(i) << rPoints[i].x, rPoints[i].y, rPoints[i].z;
	}

	std::cout << "����ת����Դ�����Ϊ��\n" << emarkerPoints << std::endl;
	std::cout << "���ڻ�׼�ľ���BΪ��\n" << esourcePoints << std::endl;



	ifstream fileDst, filedst2;
	vector<glm::dvec3> dst, dstr;
	Eigen::MatrixXd ePoints(77, 3);
	Eigen::MatrixXd mPoints(77, 3);

	int i, j;
	filedst2.open("��20ƥ��.txt");
	assert(filedst2.is_open());   //��ʧ��,�����������Ϣ,����ֹ�������� 
	while (!filedst2.eof())
	{
		for (i = 0; i <77; i++)
		{
			for (j = 0; j < 3; j++) {
				filedst2 >> mPoints(i, j);
			}
		}
	}

	std::cout << "ʵ��ת���㣺\n" << mPoints << std::endl;

	fileDst.open("��׼20.txt");   //���ļ����������ļ��������� 
	assert(fileDst.is_open());   //��ʧ��,�����������Ϣ,����ֹ�������� 
	while (!fileDst.eof())
	{
		for (i = 0; i <77; i++)
		{
			for (j = 0; j < 3; j++) {
				fileDst >> ePoints(i, j);
			}
		}
	}

	std::cout << "ʵ��ƥ���Դ��\n" << ePoints << std::endl;


	Eigen::Matrix3d R;
	Eigen::Vector3d t;
	svdRegistration(mPoints, ePoints, R, t);

	std::cout << "Rotation Matrix: \n" << R << std::endl;
	std::cout << "Translation Vector: \n" << t << std::endl;
	// �õ�rt
	
	return 0;

}




