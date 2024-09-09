#pragma once
#define _CRT_SECURE_NO_WARNINGS 1 
#include "4.0Function.h"
#include "calibFunction.h"
#include<fstream >
#include<iostream>   //定义输入输出
#include <Eigen/Dense>
#define M_PI 3.141592653589793238462643383279502884L  //定义Π

using namespace std;
using namespace cv;

//输入相机标定参数     √
void getLens(vector<SLSLens>& lens)   //填充镜头参数
{
	lens.resize(2);

	////可靠性测头0221-1714-11张图标定
	//lens[0].fc = glm::dvec2(4773.2795353944921, 4772.6264228375512);
	//lens[0].cc = glm::dvec2(1231.4497029294232, 1032.3873472430453);
	//lens[0].kc = glm::dvec4(-0.12701305511156191, 0.15541318698471046, -1.6227446198730815e-05, 0.00032514011602197769);
	//lens[1].fc = glm::dvec2(4764.5770328151948, 4763.610487964007);
	//lens[1].cc = glm::dvec2(1212.965906179596, 1004.9854018185205);
	//lens[1].kc = glm::dvec4(-0.13252868171275628,0.1840811877094049,- 0.00043647012367697741,1.2810210037806573e-06);
	//lens[1].t = glm::dvec3(-246.10274296492054, 0.28367346095095974, 34.177335703461331);
	//lens[1].om = glm::dvec3(-0.001271512517305174, 0.37359810297905127, 0.0010779583218951921);

	////可靠性测头0221-1714-9张图标定
	//lens[0].fc = glm::dvec2(4772.230255757343, 4771.8469376762732);
	//lens[0].cc = glm::dvec2(1231.5190856819736, 1032.2014852840907);
	//lens[0].kc = glm::dvec4(-0.12695758865649276, 0.14799941511556794, -9.160216377875432e-05, 0.00022453923910428439);
	//lens[1].fc = glm::dvec2(4766.1610303247335, 4765.3422571861483);
	//lens[1].cc = glm::dvec2(1213.0655522935087, 1004.6524196687323);
	//lens[1].kc = glm::dvec4(-0.13236990511535995, 0.17913426755421982, -0.00046439352147696682, -1.0503040628457219e-05);
	//lens[1].t = glm::dvec3(-246.05354404822063, 0.28782363564970914, 34.511656989066353);
	//lens[1].om = glm::dvec3(-0.0012964315619781292, 0.37358841420675437, 0.0010564013365778814);

	////可靠性测头0221-1530-11张图标定
	//lens[0].fc = glm::dvec2(4773.6189386907336, 4773.0338824235114);
	//lens[0].cc = glm::dvec2(1233.7405452279036, 1033.2225983538806);
	//lens[0].kc = glm::dvec4(-0.12772523326929291, 0.1515861851142708, -4.6761192974565691e-05, 0.0003279349328924619);
	//lens[1].fc = glm::dvec2(4767.1663345991847, 4766.143554531267);
	//lens[1].cc = glm::dvec2(1212.3156054012388, 1007.507313143862);
	//lens[1].kc = glm::dvec4(-0.13389506520082187, 0.1986125650292235, -0.0004451613549188851, -2.6267153113192434e-05);
	//lens[1].t = glm::dvec3(-246.47228943570559, 0.16429235137992493, 33.289910050813575);
	//lens[1].om = glm::dvec3(-0.00073013443714141997, 0.36961928295029728, 0.0013377298265264271);

	////可靠性测头0221-1530-9张图标定
	//lens[0].fc = glm::dvec2(4773.3215368797491, 4772.9566636006157);
	//lens[0].cc = glm::dvec2(1232.4480179542807, 1032.8404926918581);
	//lens[0].kc = glm::dvec4(-0.12717460630401003, 0.14392375552552497, -6.7751720126880567e-05, 0.00023770913786245364);
	//lens[1].fc = glm::dvec2(4766.9370258030913, 4766.2530405379584);
	//lens[1].cc = glm::dvec2(1211.7037019879865, 1006.8226078220326);
	//lens[1].kc = glm::dvec4(-0.1327340268118167, 0.18213649231635645, -0.00046330308110660904, -7.7810975671610156e-05);
	//lens[1].t = glm::dvec3(-246.45863228300553, 0.15359992460411576, 33.346397338981674);
	//lens[1].om = glm::dvec3(-0.00081970355598332001, 0.36947715532385517, 0.0013011315496112782);

	//Kasyapa0112
	lens[0].fc = glm::dvec2(4760.7340407169004,	4760.8033347408145);
	lens[0].cc = glm::dvec2(1227.3429127474881,	1005.849254096694);
	lens[0].kc = glm::dvec4(-0.13401064140210464,0.21118889464166818,- 7.6666596182125709e-05,- 3.7849380594359916e-06);
	lens[1].fc = glm::dvec2(4776.0656094764772,4775.4302990769484);
	lens[1].cc = glm::dvec2(1235.6467678692861,1019.8631712414935);
	lens[1].kc = glm::dvec4(-0.13154709038851434,0.16921970592360405,- 0.00033313902352649801,- 0.00012566200378862836);
	////标定前
	//lens[1].t = glm::dvec3(4.73735879968609, -1504.74202960518, 390.677621811001);
	//lens[1].om = glm::dvec3(-0.537923455094985, 0.00481218438872776, 0.00548677506596134);
	//标定后
	//lens[1].t = glm::dvec3(-240.2517735504145,-1.0630372350726096,45.037308244900764);
	//lens[1].om = glm::dvec3(0.001406449955126822,0.41549183274000862,0.0056949008065485269);
	lens[1].t = glm::dvec3(-239.38397468393467,- 0.97811582492819016,46.937867460667903	);
	lens[1].om = glm::dvec3(0.0011634137297956689, 0.42446252153767372,0.0050883379834097252);
}

//计算本质矩阵     √
void getEandK(vector<SLSLens>& lens, glm::dmat3& essentialMatrix, vector<glm::dmat3>& K_inv)          //计算本质矩阵和 矫正矩阵的逆矩阵
{
	//计算R
	Mat om_src = (Mat_<double>(3, 1) << lens[1].om.x, lens[1].om.y, lens[1].om.z);
	Mat R_dst;
	Rodrigues(om_src, R_dst);
	//cout << R_dst << endl;
	lens[1].R = glm::dmat3(R_dst.at<double>(0, 0), R_dst.at<double>(1, 0), R_dst.at<double>(2, 0),
		R_dst.at<double>(0, 1), R_dst.at<double>(1, 1), R_dst.at<double>(2, 1),
		R_dst.at<double>(0, 2), R_dst.at<double>(1, 2), R_dst.at<double>(2, 2));

	//cout << "lens[1].om: " << endl;
	//cout << lens[1].om.x << " " << lens[1].om.y << " " << lens[1].om.z << endl;
	//cout << "lens[1].R: " << endl;
	//cout << lens[1].R[0][0] << " " << lens[1].R[0][1] << " " << lens[1].R[0][2] << endl;
	//cout << lens[1].R[1][0] << " " << lens[1].R[1][1] << " " << lens[1].R[1][2] << endl;
	//cout << lens[1].R[2][0] << " " << lens[1].R[2][1] << " " << lens[1].R[2][2] << endl;

	//计算反对称矩阵
	glm::dmat3 antisymmetricMatrix;
	//glm::dmat3 essentialMatrix;
	antisymmetricMatrix[0][0] = 0.0f;
	antisymmetricMatrix[1][0] = -lens[1].t[2];
	antisymmetricMatrix[2][0] = lens[1].t[1];
	antisymmetricMatrix[0][1] = lens[1].t[2];
	antisymmetricMatrix[1][1] = 0.0f;
	antisymmetricMatrix[2][1] = -lens[1].t[0];
	antisymmetricMatrix[0][2] = -lens[1].t[1];
	antisymmetricMatrix[1][2] = lens[1].t[0];
	antisymmetricMatrix[2][2] = 0.0f;

	//cout << "antisymmetricMatrix: " << endl;
	//cout << antisymmetricMatrix[0][0] << "\t" << antisymmetricMatrix[0][1] << "\t" << antisymmetricMatrix[0][2] << endl;
	//cout << antisymmetricMatrix[1][0] << "\t" << antisymmetricMatrix[1][1] << "\t" << antisymmetricMatrix[1][2] << endl;
	//cout << antisymmetricMatrix[2][0] << "\t" << antisymmetricMatrix[2][1] << "\t" << antisymmetricMatrix[2][2] << endl;

	//计算本质矩阵E
	essentialMatrix = antisymmetricMatrix * lens[1].R;      //等于反对称矩阵*镜头旋转矩阵
	//cout << "\nessentialMatrix: " << endl;
	//cout << essentialMatrix[0][0] << "\t" << essentialMatrix[0][1] << "\t" << essentialMatrix[0][2] << endl;
	//cout << essentialMatrix[1][0] << "\t" << essentialMatrix[1][1] << "\t" << essentialMatrix[1][2] << endl;
	//cout << essentialMatrix[2][0] << "\t" << essentialMatrix[2][1] << "\t" << essentialMatrix[2][2] << endl;

	  //求两个相机的K_inv
	vector<Eigen::Matrix3d> KK;    //定义了一个可以储存多个对象的变量KK
	KK.resize(2);         //调整向量KK的大小是2
	for (int i = 0; i < 2; i++)
	{
		KK[i] << lens[i].fc.x, 0.0f, lens[i].cc.x, 0.0f, lens[i].fc.y, lens[i].cc.y, 0.0f, 0.0f, 1.0f;
	}  //迭代两次   将上述变量赋值到KK[i]矩阵的不同位置，将KK向量中每个Eigen::Matrix3d类型变量初始化一个值

	vector<Eigen::Matrix3d> KK_inv;  //储存KK向量中每个矩阵的逆矩阵
	KK_inv.resize(2);
	for (int i = 0; i < 2; i++)
	{
		KK_inv[i] = KK[i].inverse();
	}

	//vector<glm::dmat3> K_inv;
	K_inv.resize(2);
	for (int i = 0; i < 2; i++)
	{
		K_inv[i][0][0] = KK_inv[i](0, 0);
		K_inv[i][0][1] = KK_inv[i](1, 0);
		K_inv[i][0][2] = KK_inv[i](2, 0);
		K_inv[i][1][0] = KK_inv[i](0, 1);
		K_inv[i][1][1] = KK_inv[i](1, 1);
		K_inv[i][1][2] = KK_inv[i](2, 1);
		K_inv[i][2][0] = KK_inv[i](0, 2);
		K_inv[i][2][1] = KK_inv[i](1, 2);
		K_inv[i][2][2] = KK_inv[i](2, 2);
	}   //KK矩阵中元素复制到K矩阵中对应位置

	//cout << "K_inv(1): " << endl;
	//cout << K_inv[1][0][0] << " " << K_inv[1][0][1] << " " << K_inv[1][0][2] << endl;
	//cout << K_inv[1][1][0] << " " << K_inv[1][1][1] << " " << K_inv[1][1][2] << endl;
	//cout << K_inv[1][2][0] << " " << K_inv[1][2][1] << " " << K_inv[1][2][2] << endl;
}

//提取二维标志点
void getMarkers(cv::Mat imageL, SLSLens lens,    
	vector<glm::dvec2>& markers, vector<glm::dvec2>& undistortedMarkersL)  //getMarkers的参数有 输入的图像，镜头参数对象，用于存储检测到的标记点，去畸变后的标记点
{
	//---------------------------------------------------------提取标志点二维坐标---------------------------------------------------------
	//提取标志点
	vector<pair<glm::dvec2, float>> markersL;//定义存储标记点的向量，由二维坐标和一个float型组成
	detectMarkers(markersL, imageL);  //检测图像imageL中的标记点，存储到markersL向量中
	cout << "\n原始标志点二维坐标：" << endl;
	for (auto j = 0; j < markersL.size(); ++j) //对于markersL向量中的每个标记点
	{
		markers.push_back(markersL[j].first);   //将二维坐标添加到向量中
		cout << setprecision(15) << markers[j].x << " " << markers[j].y << endl;//输出精度为15位小数的二维坐标
	}

	//归一化
	vector<pair<glm::dvec2, double>> normalizedMarkersL;
	normalizedMarkersL.resize(markersL.size());//同一向量大小
	//vector<glm::dvec2> undistortedMarkersL;
	
	cout << "\n去畸变标志点：" << endl;
	for (auto it = 0; it < markersL.size(); ++it)
	{
		normalizedMarkersL[it].first = normalize(markersL[it].first, lens);
		undistortedMarkersL.push_back(normalizedMarkersL[it].first * lens.fc + lens.cc);
		cout << setprecision(8) << undistortedMarkersL[it].x << " " << undistortedMarkersL[it].y << endl;
	}
	//将归一化的结果存储在undistortedMarkersL向量中
}


//重建三维点云
vector<glm::dvec3> reconstruct3dPoints(vector<SLSLens> lens, cv::Mat imageL, cv::Mat imageR,
	vector<glm::dvec2>& matchMarkersL, vector<glm::dvec2>& matchMarkersR,	//左右侧图像匹配标记点
	vector<double>& epipL,vector<double>& epipR)   //引用参数
{
	matchMarkersL.clear();
	matchMarkersR.clear();
	epipL.clear();
	epipR.clear();       
	glm::dmat3 essentialMatrix;
	vector<glm::dmat3> K_inv;
	getEandK(lens, essentialMatrix, K_inv);  //调用 getEandK 函数来获取镜头参数 lens 对应的本质矩阵和逆相机矩阵

	vector<glm::dvec2> markersL, markersR, undistortedMarkersL, undistortedMarkersR;
	getMarkers(imageL, lens[0], markersL, undistortedMarkersL);
	getMarkers(imageR, lens[1], markersR, undistortedMarkersR);   //getMarkers 函数来获取图像 imageL 和 imageR 中的标记点

	//---------------------------------------------------------重建三维点---------------------------------------------------------
	vector<glm::dvec3> markerPoints;
	//vector<glm::dvec2> matchedMarkersL;
	//vector<glm::dvec2> matchedMarkersR;

	//输出
	std::ofstream File(".\\markerEpipolarError.txt", std::ios::ate);//创建输出文件
	File.precision(8);//八位有效数字
	File << "This is a new set of markers." << std::endl;

	double minDepth = 1.0;
	double maxDepth = 100000.0;//5000.0
	cout << "----------------------------重建三维点----------------------------" << endl;

	//遍历左相机第一张图的每个点
	for (auto i = 0; i < markersL.size(); ++i)
	{
		//cout << "左：" << markersL[0][i].first.x << " " << markersL[0][i].first.y << endl;
		vector<glm::dvec2> targetMarkers;                  //存储点的二维坐标
		vector<double> targetEpipLs,targetEpipRs;           //存储图像上的极线方程

		//遍历右相机第一张图的每个点
		cout << "\n点" << i << endl;
		cout << "左：" << markersL[i].x <<" " << markersL[i].y << endl;
		for (auto j = 0; j < markersR.size(); ++j)
		{
			//计算极线误差和匹配误差
			glm::dvec2 tem = epipolarErr(undistortedMarkersL[i], undistortedMarkersR[j], essentialMatrix, K_inv);
			glm::dvec3 epipolarLine;
			glm::dvec2 tem2 = epipolarErr2(undistortedMarkersL[i], undistortedMarkersR[j], epipolarLine, essentialMatrix, K_inv);//2023/08/16 WYX
			
			double fileOutputEpipolarErrThreshold = 1.0;            //0.25f;//mzy20231027 输出极线匹配误差的阈值
			if (tem.x / tem.y < fileOutputEpipolarErrThreshold)
			{
				File << "imgL(" << markersL[i].x << ", " << markersL[i].y
					<< ") imgR(" << markersR[j].x << ", " << markersR[j].y << ") ";
				File /*<< normalizedMarkers[0][i].first.x << " " << normalizedMarkers[0][i].first.y << " "
					<< normalizedMarkers[1][j].first.x << " " << normalizedMarkers[1][j].first.y << " "*/
					<< "   epipolarErrL2R = " << (tem.x / tem.y)
					<< "   epipolarErrR2L = " << (tem2.x / tem2.y) << std::endl;
			}


			//保留极线匹配误差小于阈值的点
			if (checkEpipolar(undistortedMarkersL[i],
				undistortedMarkersR[j], fileOutputEpipolarErrThreshold, essentialMatrix, K_inv))
			{
				targetMarkers.push_back(markersR[j]);
				targetEpipLs.push_back(tem.x / tem.y);
				targetEpipRs.push_back(tem2.x / tem2.y);
			}
		}

		//确定匹配点
		glm::dvec2 matchMarker;    //存储匹配点的二维坐标
		double matchEpipL, matchEpipR;   //存储最佳
		matchEpipL = 100.0;
		if (targetMarkers.size() != 0)
		{
			double corrCoeff;   //存储相关系数
			glm::dvec3 point = glm::dvec3();
			int roiPara = 50;
			if (targetMarkers.size() > 1)
			{
				for (int markers_i = 0; markers_i < targetMarkers.size(); ++markers_i)
				{
					if (targetEpipLs[markers_i] < matchEpipL)  //判断效果是否更佳
					{
						matchMarker = targetMarkers[markers_i];//更新匹配标记点
						matchEpipL = targetEpipLs[markers_i];
						matchEpipR = targetEpipRs[markers_i];  //更新匹配点在图像上的极限方程参数
					}
				}
				////根据灰度值寻找唯一匹配点
				//corrCoeff = checkMarkers(imageL, imageR, markersL[i], targetMarkers, targetEpipLs, targetEpipRs,
				//	matchMarker, matchEpipL, matchEpipR, roiPara);
			}

			else
			{
				matchMarker = targetMarkers[0];
				matchEpipL = targetEpipLs[0];
				matchEpipR = targetEpipRs[0];
			}
				
			cout << "右：" << matchMarker.x << " " << matchMarker.y << endl << endl;

			//point = reconstructor.triangulate2(
			//	glm::dvec2(markers[0][i].first), glm::dvec2(matchMarker.first), true); //2023/08/18 WYX 如果用TopCam，用triangulate2重建

			//三角重建
			//cout << "------------------------" << i << "----------------------------" << endl;
			//cout << "markersL[i]: " << markersL[i].x << " " << markersL[i].y << endl;
			//cout << "matchMarker: " << matchMarker.x << " " << matchMarker.y << endl;
			point = triangulate(glm::dvec2(markersL[i]), matchMarker, lens);//2023/08/18 WYX 如果用竖柱子，用triangulate重建
			cout << "三维点："<<point.x << " " << point.y << " " << point.z << endl;

			if (point.z > minDepth && point.z < maxDepth)
			{
				//应该和绘制加锁
				//markerPoints.push_back(glm::dvec3(reconstructor.getGlobalPose() * glm::dvec4(point, 1.0f)));				//2022/11/05 YCC：乘上 globalpose
				markerPoints.push_back(point);
				matchMarkersL.push_back(glm::dvec2(markersL[i]));
				matchMarkersR.push_back(glm::dvec2(matchMarker));
				epipL.push_back(matchEpipL);
				epipR.push_back(matchEpipR);
			}
		}
	}

	File << std::endl;
	File.close();

	//cout << "\n重建三维点（未偏心校正）：" << endl;
	//for (int it = 0; it < markerPoints.size(); it++)
	//{
	//	cout << markerPoints[it].x << " " << markerPoints[it].y << " " << markerPoints[it].z << endl;
	//}

	return markerPoints;
}

//得到双目相机系统总体姿态
string calibRT(vector<glm::dvec3> source,	vector<glm::dvec3> target,      //源点的三维坐标，目标点的三维坐标
	const glm::dmat4& rt, 	vector<glm::dvec3> pt3d,                        //旋转平移矩阵，计算得到的三维点坐标
	vector<glm::dvec2> imgL, vector<glm::dvec2> imgR,                        //图像中二维点坐标
	vector<double> epipL, vector<double> epipR,                      //图像中极线方程参数
	vector<SLSLens>& lens, bool doCalib)                        //bool变量用于指示是否进行标定
{
	std::vector<int> index;     // 选出与激光跟踪仪坐标点匹配的标志点
	for (int i = 0; i < source.size(); i++)
	{
		for (int j = 0; j < pt3d.size(); j++)
		{
			if (glm::length(source[i] - pt3d[j]) < 2)   //判断点源与计算得到的三维坐标的距离值
				index.push_back(j);
		}
	}
	//OutputDebugStringA(("index size = " + std::to_string(index.size()) + "\r\n").c_str());
	if (index.size() < 6)
	{
		return "标志点少于6个！";
	}
	
	cv::Mat rmat, tvec, rvec;

	if (doCalib)
	{
		//Solvpnp计算外参数
		std::vector<cv::Point2d> imagePointsLeft;     // 获取世界坐标系三维坐标、左右图像二维坐标
		std::vector<cv::Point2d> imagePointsRight;
		std::vector<cv::Point3d> objectPoints;
		for (int i = 0; i < index.size(); i++)
		{
			imagePointsLeft.emplace_back(imgL[index[i]].x, imgL[index[i]].y);
			imagePointsRight.emplace_back(imgR[index[i]].x, imgR[index[i]].y);  //获取左右图像中点的xy坐标，添加到imagepoint.函数中
			objectPoints.emplace_back(target[i].x, target[i].y, target[i].z);  //从target中获取点的三维坐标
		}

		cv::Mat cameraMatrixLeft = (cv::Mat_<double>(3, 3) << lens[0].fc.x, 0, lens[0].cc.x,
			0, lens[0].fc.y, lens[0].cc.y,0, 0, 1);                      //用lens[0]中的数值初始化相机矩阵
		cv::Mat cameraMatrixRight = (cv::Mat_<double>(3, 3) << lens[1].fc.x, 0, lens[1].cc.x,
			0, lens[1].fc.y, lens[1].cc.y,0, 0, 1);
		cv::Mat distCoeffsLeft = (cv::Mat_<double>(4, 1) << lens[0].kc.x, lens[0].kc.y, lens[0].kc.z, lens[0].kc.w);   //畸变系数
		cv::Mat distCoeffsRight = (cv::Mat_<double>(4, 1) << lens[1].kc.x, lens[1].kc.y, lens[1].kc.z, lens[1].kc.w);
		cv::Mat rvecLeft, rvecRight, tvecLeft, tvecRight;            //存储旋转、平移向量
		cv::Mat rvecLeft1, rvecRight1, tvecLeft1, tvecRight1;          
		cv::Mat inliers0, inliers1;                  //存储内点


		int test = 0;
		if (test == 0)
		{
			//相机图像点imagepoint，相机矩阵，畸变系数
			cv::solvePnP(objectPoints, imagePointsLeft, cameraMatrixLeft, distCoeffsLeft, rvecLeft, tvecLeft, false, cv::SOLVEPNP_ITERATIVE);
			cv::solvePnP(objectPoints, imagePointsRight, cameraMatrixRight, distCoeffsRight, rvecRight, tvecRight, false, cv::SOLVEPNP_ITERATIVE);
			//再次使用cv::solvePnP函数，对初始值进行迭代优化
			cv::solvePnP(objectPoints, imagePointsLeft, cameraMatrixLeft, distCoeffsLeft, rvecLeft, tvecLeft, true, cv::SOLVEPNP_ITERATIVE);
			cv::solvePnP(objectPoints, imagePointsRight, cameraMatrixRight, distCoeffsRight, rvecRight, tvecRight, true, cv::SOLVEPNP_ITERATIVE);
		}
		else if (test == 1)  // 先EPNP计算初值，再进行Interative
		{
			//图像点，物体点，相机矩阵，畸变系数，采用EPnP算法
			cv::solvePnP(objectPoints, imagePointsLeft, cameraMatrixLeft, distCoeffsLeft, rvecLeft, tvecLeft, false, cv::SOLVEPNP_EPNP);
			cv::solvePnP(objectPoints, imagePointsRight, cameraMatrixRight, distCoeffsRight, rvecRight, tvecRight, false, cv::SOLVEPNP_EPNP);
			//再次使用cv::，迭代优化
			cv::solvePnP(objectPoints, imagePointsLeft, cameraMatrixLeft, distCoeffsLeft, rvecLeft, tvecLeft, true, cv::SOLVEPNP_ITERATIVE);
			cv::solvePnP(objectPoints, imagePointsRight, cameraMatrixRight, distCoeffsRight, rvecRight, tvecRight, true, cv::SOLVEPNP_ITERATIVE);
		}
		if (test == 0) // RANSAC SOLVEPNP
		{
			//SOLVEPNPRANSAC阈值
			//bool useExtrinsicGuess = false, int iterationsCount = 100, loat reprojectionError = 8.0, double confidence = 0.99, 
			//OutputArray inliers = noArray(), int flags = SOLVEPNP_ITERATIVE
			cv::solvePnPRansac(objectPoints, imagePointsLeft, cameraMatrixLeft, distCoeffsLeft, rvecLeft1, tvecLeft1, false,
				100000, 0.1, 0.999, inliers0, cv::SOLVEPNP_ITERATIVE);
			cv::solvePnPRansac(objectPoints, imagePointsRight, cameraMatrixRight, distCoeffsRight, rvecRight1, tvecRight1, false,
				100000, 0.1, 0.999, inliers1, cv::SOLVEPNP_ITERATIVE);
		}

		cv::Mat rmatLeft, rmatRight;
		cv::Rodrigues(rvecLeft, rmatLeft);
		cv::Rodrigues(rvecRight, rmatRight);
		//cv::Mat rmat, tvec, rvec;
		rmat = rmatRight * rmatLeft.t();
		tvec = -rmatRight * rmatLeft.t() * tvecLeft + tvecRight;
		cv::Rodrigues(rmat, rvec);

		//调试SolvePnPRansac
		cv::Mat rmatLeft1, rmatRight1;
		cv::Rodrigues(rvecLeft1, rmatLeft1);
		cv::Rodrigues(rvecRight1, rmatRight1);
		cv::Mat rmat1, tvec1, rvec1;
		rmat1 = rmatRight1 * rmatLeft1.t();
		tvec1 = -rmatRight1 * rmatLeft1.t() * tvecLeft1 + tvecRight1;
		cv::Rodrigues(rmat1, rvec1);
		//矩阵运算，得到双目相机系统总体姿态
	}

	std::string tmpName = "SolvePnPResult.txt";
	if(!doCalib) tmpName= "SolvePnPResult验证.txt";
	//std::string tmpName = "./solvePnP/OutFiles/SolvePnPResult.txt";
	std::ofstream oFile(tmpName, std::ios::ate);
	oFile.precision(15);           //15位小数

	oFile << "当前相机外参数RT：" << "\n"
		<< "    om: [" << std::setprecision(15) << lens[1].om.x << ", "
		<< std::setprecision(15) << lens[1].om.y << ", "
		<< std::setprecision(15) << lens[1].om.z << "]\n"
		<< "    t: [" << std::setprecision(15) << lens[1].t.x << ", "
		<< std::setprecision(15) << lens[1].t.y << ", "
		<< std::setprecision(15) << lens[1].t.z << "]\n";

	//oFile << "计算的双目RTCamtoGlobal（SolvePnP）：" << "\n"
	//	<< "    omR: [" << std::setprecision(15) << rvecRight.ptr<double>(0)[0] << ", "
	//	<< std::setprecision(15) << rvecRight.ptr<double>(1)[0] << ", " 
	//	<< std::setprecision(15) << rvecRight.ptr<double>(2)[0] << "]\n"
	//	<< "    tR: [" << std::setprecision(15) << tvecRight.ptr<double>(0)[0] << ", "
	//	<< std::setprecision(15) << tvecRight.ptr<double>(1)[0] << ", "
	//	<< std::setprecision(15) << tvecRight.ptr<double>(2)[0] << "]\n"

	//	<< "    omL: [" << std::setprecision(15) << rvecLeft.ptr<double>(0)[0] << ", "
	//	<< std::setprecision(15) << rvecLeft.ptr<double>(1)[0] << ", "
	//	<< std::setprecision(15) << rvecLeft.ptr<double>(2)[0] << "]\n"
	//	<< "    tL: [" << std::setprecision(15) << tvecLeft.ptr<double>(0)[0] << ", "
	//	<< std::setprecision(15) << tvecLeft.ptr<double>(1)[0] << ", "
	//	<< std::setprecision(15) << tvecLeft.ptr<double>(2)[0] << "]\n";

	//oFile << "\n-------------------------------------------------------- "
	//	<< "\n计算的双目外参数RT（SolvePnPRansac）：" << "\n"
	//	<< "    om: [" << std::setprecision(15) << rvec1.ptr<double>(0)[0] << ", "
	//	<< std::setprecision(15) << rvec1.ptr<double>(1)[0] << ", "
	//	<< std::setprecision(15) << rvec1.ptr<double>(2)[0] << "]\n"
	//	<< "    t: [" << std::setprecision(15) << tvec1.ptr<double>(0)[0] << ", "
	//	<< std::setprecision(15) << tvec1.ptr<double>(1)[0] << ", "
	//	<< std::setprecision(15) << tvec1.ptr<double>(2)[0] << "]\n";

	//oFile << "\nRansac： 左相机内点 num:" << inliers0.rows << "\n";
	//for (int i = 0; i < inliers0.rows; i++)
	//{
	//	for (int j = 0; j < inliers0.cols; j++)
	//	{
	//		oFile << std::setprecision(15) << inliers0.ptr<int>(i)[j] + 1 << " ";
	//	}
	//}
	//oFile << "\nRansac： 右相机内点 num:" << inliers1.rows << "\n";
	//for (int i = 0; i < inliers1.rows; i++)
	//{
	//	for (int j = 0; j < inliers1.cols; j++)
	//	{
	//		oFile << std::setprecision(15) << inliers1.ptr<int>(i)[j] + 1 << " ";
	//	}
	//}
	oFile << "\n二维点坐标： Left Camera --  Right Camera " << "\n";
	for (int i = 0; i < index.size(); i++)
	{
		oFile << " [" << i + 1 << "]: " << std::setprecision(15) << imgL[index[i]].x << " "
			<< std::setprecision(15) << imgL[index[i]].y << "  "
			<< std::setprecision(15) << imgR[index[i]].x << " "
			<< std::setprecision(15) << imgR[index[i]].y << "\n";
	}
	//将imgL和imgR数组中的二维坐标写入文件，每行表示一个点及其左右相机坐标

	oFile << "\n极线匹配误差： Left2Right  --  Right2Left " << "\n";
	for (int i = 0; i < index.size(); i++)
	{
		oFile << " [" << i + 1 << "]: "
			<< std::setprecision(15) << epipL[index[i]] << " "
			<< std::setprecision(15) << epipR[index[i]] << "\n";
	}

	oFile << "\n重建三维点坐标： Measure Markers Point3d " << "\n";
	for (int i = 0; i < index.size(); i++)
	{
		oFile << " [" << i + 1 << "]: " << std::setprecision(15) << source[i].x << " "
			<< std::setprecision(15) << source[i].y << " "
			<< std::setprecision(15) << source[i].z << "\n";
	}

	oFile << "\n激光跟踪仪三维点坐标： Laser Markers Point3d " << "\n";
	for (int i = 0; i < index.size(); i++)
	{
		oFile << " [" << i + 1 << "]: " << std::setprecision(15) << target[i].x << " "
			<< std::setprecision(15) << target[i].y << " "
			<< std::setprecision(15) << target[i].z << "\n";
	}

	oFile << "\n三维点匹配误差： MatchAfter Error " << "\n";
	for (int i = 0; i < index.size(); i++)
	{
		glm::dvec3 tmp = glm::dvec3(rt * glm::dvec4(target[i], 1));
		oFile << " [" << i + 1 << "]: " << std::setprecision(15) << glm::length(source[i] - tmp)
			<< "    " << std::setprecision(15) << source[i].x - tmp.x << " "
			<< std::setprecision(15) << source[i].y - tmp.y << " "
			<< std::setprecision(15) << source[i].z - tmp.z << " " << "\n";
	}

	if (doCalib)
	{
		oFile << "\n计算的双目外参数RT（SolvePnP）：" << "\n"
			<< "    om: [" << std::setprecision(15) << rvec.ptr<double>(0)[0] << ", "
			<< std::setprecision(15) << rvec.ptr<double>(1)[0] << ", "
			<< std::setprecision(15) << rvec.ptr<double>(2)[0] << "]\n"
			<< "    t: [" << std::setprecision(15) << tvec.ptr<double>(0)[0] << ", "
			<< std::setprecision(15) << tvec.ptr<double>(1)[0] << ", "
			<< std::setprecision(15) << tvec.ptr<double>(2)[0] << "]\n";

		lens[1].om.x = rvec.ptr<double>(0)[0];
		lens[1].om.y = rvec.ptr<double>(1)[0];
		lens[1].om.z = rvec.ptr<double>(2)[0];
		lens[1].t.x = tvec.ptr<double>(0)[0];
		lens[1].t.y = tvec.ptr<double>(1)[0];
		lens[1].t.z = tvec.ptr<double>(2)[0];
	}

	oFile.close();

	/*std::string matchedMarkersFileName = "CalibmatchedMarkers.txt";    // Solvepnp需要的三维坐标，二维坐标
	std::ofstream outputFile(matchedMarkersFileName, std::ios::ate);
	for (int i = 0; i < index.size(); i++)
	{
		std::setprecision(15);
		outputFile << std::setprecision(15)
				   << target[i].x << " " << target[i].y << " " << target[i].z << " "
				   << imgL[index[i]].x << " " << imgL[index[i]].y << " "
				   << imgR[index[i]].x << " " << imgR[index[i]].y << "\n";
	}	*/

	return std::string();
}


//双精度浮点数4*4矩阵
glm::dmat4 matchPointSet2(const vector<glm::dvec3>& sourcePoints, 	vector<glm::dvec3>& targetPoints, //源点集合就是待匹配点集合
	vector<glm::dvec3>& matchedSourcePoints, vector<glm::dvec3>& matchedTargetPoints)  //成功匹配的源点集合
{
	//OutputDebugStringA(std::string("sourcePoints" + std::to_string(sourcePoints.size()) + "\r\n").c_str());
	//OutputDebugStringA(std::string("targetPoints" + std::to_string(targetPoints.size()) + "\r\n").c_str());
	if (sourcePoints.size() < 3 || targetPoints.size() < 3)
		return glm::dmat4();

	matchedSourcePoints.clear();
	matchedTargetPoints.clear();//清空集合里的元素

	//SLSClock clock;
	//clock.begin();
	std::vector<std::tuple<int, int, int, double>> sourceDistances;              
	sourceDistances.clear();
	std::vector<std::tuple<int, int, int, double>> targetDistances;              //清空向量元素
	targetDistances.clear();
	for (auto i = 0; i < sourcePoints.size(); ++i)
		for (auto j = i + 1; j < sourcePoints.size(); ++j)
			sourceDistances.push_back(std::make_tuple(0, i, j,
				glm::distance(sourcePoints[i], sourcePoints[j])));          //计算每对点之间的距离  0表示源点集合，1表示目标点集合
	for (auto i = 0; i < targetPoints.size(); ++i)
		for (auto j = i + 1; j < targetPoints.size(); ++j)
			targetDistances.push_back(std::make_tuple(1, i, j,
				glm::distance(targetPoints[i], targetPoints[j])));             //make_tuple表示包含ij和距离索引的元组
	
	auto sortDistance = [](const std::tuple<int, int, int, double>& first,
		const std::tuple<int, int, int, double>& second)
	{
		return std::get<3>(first) < std::get<3>(second);      //根据距离值对元组进行排序
	};

	std::vector<std::tuple<int, int, int, double>> distances(sourceDistances.size() +
		targetDistances.size());        //存储合并后的距离 元组
	std::merge(sourceDistances.begin(), sourceDistances.end(), targetDistances.begin(),
		targetDistances.end(), distances.begin());      //合并
	std::sort(distances.begin(), distances.end(), sortDistance); //排序
	std::vector<double> distanceCenters;   //存储距离的中心位置
	distanceCenters.clear();
	auto currentDistance = 0.0f;
	double matchPointSetThreshold = 15.0;
	
	for (const auto& distance : distances)
		if (fabs(std::get<3>(distance) - currentDistance) > matchPointSetThreshold)
		{
			//判断第四个元素的距离值是否大于15
			distanceCenters.push_back(std::get<3>(distance));
			currentDistance = std::get<3>(distance);  //更新
		}
	//筛选出一组距离值
	
	//printf("%d %d\n", distances.size(), distanceCenters.size());

	//clock.begin();
	std::vector<std::vector<std::tuple<int, int, int, double>>> distancesList(distanceCenters.size());
	//将距离元组按照中心距离进行分组
	for (auto i = 0; i < distanceCenters.size(); ++i)
		for (const auto& distance : distances)
			if (fabs(std::get<3>(distance) -
				distanceCenters[i]) < matchPointSetThreshold)
				distancesList[i].push_back(distance);  //将该距离元组添加到向量中

	for (auto& distances : distancesList)//循环迭代向量中的每组距离元组
	{
		if (distances.size() == 1)   //判断组中是否只有一个距离元组
			distances.clear();
		else
		{
			auto findSource = false;
			auto findTarget = false;
			for (const auto& distance : distances)  //循环当前组内的距离元组
			{
				if (std::get<0>(distance) == 0)
					findSource = true;
				if (std::get<0>(distance) == 1)
					findTarget = true;
			}
			//检查每一个元素
			if (!findSource || !findTarget)
				distances.clear();
		}
	}

	//printf("%d\n", distancesList.size());

	//clock.begin();将每个距离tuple与distancesIndices向量中相应的索引相关联
	std::vector<std::vector<glm::ivec2>> distancesIndices(glm::max(
		sourcePoints.size(), targetPoints.size()));//2D向量
	for (auto i = 0; i < distancesList.size(); ++i)
	{
		for (const auto& distance : distancesList[i])
		{
			distancesIndices[std::get<1>(distance)].push_back(glm::ivec2(
				std::get<0>(distance), i));
			distancesIndices[std::get<2>(distance)].push_back(glm::ivec2(
				std::get<0>(distance), i));
		}//根据特定条件存储距离tuple的index
	}

	//clock.begin();
	std::vector<std::vector<int>> sourceDistancesIndices(sourcePoints.size());
	std::vector<std::vector<int>> targetDistancesIndices(targetPoints.size());
	for (auto i = 0; i < distancesIndices.size(); ++i)
		if (distancesIndices[i].size() >= 2)
			for (const auto& distanceIndex : distancesIndices[i])
			{
				if (distanceIndex.x == 0)
					sourceDistancesIndices[i].push_back(distanceIndex.y);
				if (distanceIndex.x == 1)
					targetDistancesIndices[i].push_back(distanceIndex.y);
			}

	for (auto& sourceDistanceIndices : sourceDistancesIndices)
		std::sort(sourceDistanceIndices.begin(), sourceDistanceIndices.end());
	for (auto& targetDistanceIndices : targetDistancesIndices)
		std::sort(targetDistanceIndices.begin(), targetDistanceIndices.end());

	auto findCommonDistanceIndices = [](const std::vector<int>& indices1,
		const std::vector<int>& indices2)
	{
		auto i = 0;
		auto j = 0;
		auto count = 0;
		while (i < indices1.size() && j < indices2.size())
		{
			if (indices1[i] < indices2[j])
				++i;
			else if (indices1[i] == indices2[j])
			{
				++count;
				++i;
				++j;
			}
			else
				++j;
		}
		return count;
	};

	std::vector<glm::ivec3> matchedIndices;
	matchedIndices.clear();
	for (auto i = 0; i < sourceDistancesIndices.size(); ++i)
	{
		//auto maxCount = 0;
		//auto maxTargetIndex = -1;
		std::vector<std::pair<int, int>> TargetIndices;
		TargetIndices.clear();
		TargetIndices.push_back(std::pair<int, int>(-1, 1));
		for (auto j = 0; j < targetDistancesIndices.size(); ++j)
		{
			int count = findCommonDistanceIndices(sourceDistancesIndices[i],
				targetDistancesIndices[j]);
			//可能有多个相等的最大值

			if (count >= TargetIndices.back().second)
			{
				TargetIndices.push_back(std::pair<int, int>(j, count));
			}
			//if (count > maxCount)
			//{
			//    maxCount = count;
			//    maxTargetIndex = j;
			//}
			//if (count > 0)
			//	printf("%d %d: %d\n", i, j, count);
		}
		if (TargetIndices.size() > 1)
		{
			for (auto TargetIndice : TargetIndices)
			{
				if (TargetIndice.first != -1)
				{
					matchedIndices.emplace_back(i, TargetIndice.first, TargetIndice.second);
				}
			}
		}
		//if (maxTargetIndex >= 0)
		//	matchedIndices.emplace_back(i, maxTargetIndex, maxCount);
		//printf("%d %d: %d\n", i, maxTargetIndex, maxCount);
	}

	//for (int i = 0; i < matchedIndices.size(); ++i)
	//{
	//    OutputDebugStringA(std::string("(" + std::to_string(matchedIndices[i].x) + "," + 
	//                                    std::to_string(matchedIndices[i].y) + "," + 
	//                                    std::to_string(matchedIndices[i].z) + ")\r\n").c_str());
	//}

	if (matchedIndices.size() < 3)
	{
		OutputDebugStringA(std::string("1061 return\r\n").c_str());
		return glm::dmat4();
	}

	auto checkInner = [&matchedIndices, &sourcePoints, &targetPoints]
	(const glm::dmat4& transformation, double& matchPointSetThreshold, const int& index)
	{
		return glm::distance(glm::dvec3(transformation *
			glm::dvec4(sourcePoints[matchedIndices[index].x], 1.0)),
			targetPoints[matchedIndices[index].y]) < matchPointSetThreshold;
	};

	std::vector<glm::dvec3> sourceMatchedPoints;
	std::vector<glm::dvec3> targetMatchedPoints;

	auto extractMatchedPoints = [&matchedIndices, &sourcePoints, &targetPoints,
		&sourceMatchedPoints, &targetMatchedPoints, &checkInner]
		(const glm::dmat4& transformation, double& matchPointSetThreshold)
	{
		sourceMatchedPoints.clear();
		targetMatchedPoints.clear();
		for (auto i = 0; i < matchedIndices.size(); ++i)
		{
			if (checkInner(transformation, matchPointSetThreshold, i))
			{
				sourceMatchedPoints.push_back(sourcePoints[matchedIndices[i].x]);
				targetMatchedPoints.push_back(targetPoints[matchedIndices[i].y]);
				//	printf("inner %d %d: %d\n", matchedIndices[i].x, matchedIndices[i].y,
				//		matchedIndices[i].z);
			}
			//else
			//	printf("outer %d %d: %d\n", matchedIndices[i].x, matchedIndices[i].y,
			//		matchedIndices[i].z);
		}
	};

	
	//可以穷举
	std::vector<glm::ivec3> sourceMatchedIndices;
	sourceMatchedIndices.clear();
	std::vector<glm::ivec3> targetMatchedIndices;
	targetMatchedIndices.clear();
	for (auto l = 0; l < matchedIndices.size(); ++l)
		for (auto m = l + 1; m < matchedIndices.size(); ++m)
			for (auto n = m + 1; n < matchedIndices.size(); ++n)
			{
				sourceMatchedIndices.emplace_back(matchedIndices[l].x,
					matchedIndices[m].x, matchedIndices[n].x);
				targetMatchedIndices.emplace_back(matchedIndices[l].y,
					matchedIndices[m].y, matchedIndices[n].y);
			}

	sourceMatchedPoints.resize(3);
	targetMatchedPoints.resize(3);

	auto threadsNum = omp_get_max_threads();
	std::vector<int> maxInnersNums(threadsNum, 0);
	std::vector<glm::dmat4> innerTransformations(threadsNum);

#pragma omp parallel for schedule (dynamic, 1)
	for (auto i = 0; i < sourceMatchedIndices.size(); ++i)
	{
		std::vector<glm::dvec3> _sourceMatchedPoints;
		_sourceMatchedPoints.clear();
		_sourceMatchedPoints.resize(3, glm::dvec3(0, 0, 0));
		std::vector<glm::dvec3> _targetMatchedPoints;
		_targetMatchedPoints.clear();
		_targetMatchedPoints.resize(3, glm::dvec3(0, 0, 0));
		for (auto j = 0; j < 3; ++j)
		{
			_sourceMatchedPoints[j] = sourcePoints[sourceMatchedIndices[i][j]];
			_targetMatchedPoints[j] = targetPoints[targetMatchedIndices[i][j]];
		}
		auto transformation = estimateRigidTransformation(_sourceMatchedPoints.size(),
			_sourceMatchedPoints, _targetMatchedPoints);
		if (i == 1)
		{
			OutputDebugStringA("source:\r\n");
			for (auto j = 0; j < 3; ++j)
			{
				OutputDebugStringA(std::string("(" + std::to_string(_sourceMatchedPoints[j][0]) + ",").c_str());
				OutputDebugStringA(std::string(std::to_string(_sourceMatchedPoints[j][1]) + ",").c_str());
				OutputDebugStringA(std::string(std::to_string(_sourceMatchedPoints[j][2]) + ")\r\n").c_str());
			}
			OutputDebugStringA("target:\r\n(");
			for (auto j = 0; j < 3; ++j)
			{
				OutputDebugStringA(std::string(std::to_string(_targetMatchedPoints[j][0]) + ",").c_str());
				OutputDebugStringA(std::string(std::to_string(_targetMatchedPoints[j][1]) + ",").c_str());
				OutputDebugStringA(std::string(std::to_string(_targetMatchedPoints[j][2]) + ")\r\n").c_str());
			}
			OutputDebugStringA(std::string(glm::to_string(transformation) + "\r\n").c_str());
		}
		auto checkInner_ = [ &matchedIndices, &sourcePoints, &targetPoints]
		(const glm::dmat4& transformation, double& matchPointSetThreshold, const int& index)
		{
			auto d = glm::distance(glm::dvec3(transformation *
				glm::dvec4(sourcePoints[matchedIndices[index].x], 1.0)),
				targetPoints[matchedIndices[index].y]);
			//OutputDebugStringA(std::string("dis=" + std::to_string(d) + "\tindex" + std::to_string(index) + "\r\n").c_str());
			return d < matchPointSetThreshold;
		};

		auto innersNum = 0;
		for (auto j = 0; j < matchedIndices.size(); ++j)
			if (checkInner_(transformation, matchPointSetThreshold,j))
				++innersNum;
		auto threadIndex = omp_get_thread_num();
		//OutputDebugStringA(std::string("threadIndex:" + std::to_string(threadIndex) + "\tnum:" + std::to_string(innersNum) + "\r\n").c_str());
		if (innersNum > maxInnersNums[threadIndex])
		{
			//OutputDebugStringA(std::string("threadIndex:" + std::to_string(threadIndex) + "\tnum:" + std::to_string(innersNum)).c_str());
			maxInnersNums[threadIndex] = innersNum;
			innerTransformations[threadIndex] = transformation;
		}
	}
	auto maxInnersNum = 0;
	glm::dmat4 innerTransformation;
	for (auto i = 0; i < threadsNum; ++i)
		if (maxInnersNums[i] > maxInnersNum)
		{
			maxInnersNum = maxInnersNums[i];
			innerTransformation = innerTransformations[i];
		}

	if (maxInnersNum < 3)
	{
		OutputDebugStringA(std::string("1205 return\r\n").c_str());
		return glm::dmat4();
	}


	extractMatchedPoints(innerTransformation, matchPointSetThreshold);
	

	auto transformation = estimateRigidTransformation(sourceMatchedPoints.size(),
		sourceMatchedPoints, targetMatchedPoints);

	//printf("matchedPointsNum: %zd\n", sourceMatchedPoints.size());

	if (sourceMatchedPoints.size() >= 3)
		transformation = estimateRigidTransformation(sourceMatchedPoints.size(),
			sourceMatchedPoints, targetMatchedPoints);
	//SLS::displayMatrix("matchPointSet", transformation);

	//clock.end();
	//clock.displayInterval("matchPointSet");

	//输出匹配误差 2023/06/08 WYX
	//get local time
	auto now = std::chrono::system_clock::now();
	std::time_t time = std::chrono::system_clock::to_time_t(now);
	std::tm* nowTime = std::localtime(&time);
	//write txt
	std::ofstream errorWriter;
	errorWriter.open("./markerMatchError.txt", std::ios::app);
	errorWriter << "[" << std::setw(2) << nowTime->tm_hour << ":"
		<< std::setw(2) << nowTime->tm_min << ":"
		<< std::setw(2) << nowTime->tm_sec << "]\n";
	errorWriter << "number\t" << "source pt\t\t\t" << "target pt\t\t\t\t" << "error" << std::endl;

	for (auto i = 0; i < sourceMatchedPoints.size(); ++i)
	{
		auto error = glm::distance(glm::dvec3(transformation *
			glm::dvec4(sourceMatchedPoints[i], 1.0)), targetMatchedPoints[i]);
		//printf("%d: %f\n", i, error);

		//2023/06/08 WYX
		errorWriter << i << "\t("
			<< std::setprecision(15) << (transformation * glm::dvec4(sourceMatchedPoints[i], 1.0))[0] << ", "
			<< std::setprecision(15) << (transformation * glm::dvec4(sourceMatchedPoints[i], 1.0))[1] << ", "
			<< std::setprecision(15) << (transformation * glm::dvec4(sourceMatchedPoints[i], 1.0))[2] << ")\t("
			<< std::setprecision(15) << targetMatchedPoints[i].x << ", "
			<< std::setprecision(15) << targetMatchedPoints[i].y << ", "
			<< std::setprecision(15) << targetMatchedPoints[i].z << ")\t\t"
			<< error << std::endl;
	}
	errorWriter << std::endl;
	errorWriter.close();

	for (int n = 0; n < targetPoints.size(); n++)
	{
		for (int m = 0; m < targetMatchedPoints.size(); m++)
		{
			if (targetMatchedPoints[m].x == targetPoints[n].x)
			{
				matchedTargetPoints.push_back(targetMatchedPoints[m]);
				matchedSourcePoints.push_back(sourceMatchedPoints[m]);
			}
		}
	}
	/*for (const auto &pt:sourceMatchedPoints)
	{
		matchedSourcePoints.push_back(pt);
	}

	for (const auto& pt : targetMatchedPoints)
	{
		matchedTargetPoints.push_back(pt);
	}*/

	//2023/06/12 WYX 再进行一次ICP
	//double currentMaxDistance = paras.maxDistance;
	//paras.maxDistance = paras.matchPointSetThreshold * 0.667;
	//bool currentCheckNormal = paras.checkNormal;
	//paras.checkNormal = false;

	//SLSGeometry srcGeometry;
	//SLSGeometry* tgtGeometry = new SLSGeometry;
	//std::vector<SLSGeometry*> tgtGeometries;
	//std::vector<glm::dvec3> srcPts, tgtPts;
	//for(const auto& pt : sourcePoints)
	//	srcGeometry.points.push_back(glm::dvec3(transformation * glm::dvec4(pt, 1.0)));
	//for (const auto& pt : targetPoints)
	//	tgtGeometry->points.push_back(pt);
	//tgtGeometries.push_back(std::move(tgtGeometry));
	//auto [transICP, overlapped] = iterativeClosestPoint(srcGeometry, tgtGeometries);
	//transformation = transICP * transformation;

	//paras.maxDistance = currentMaxDistance;
	//paras.checkNormal = currentCheckNormal;

	////==================== debug ====================
	////再输出一次ICP之后的偏差，检验ICP是否有用
	//now = std::chrono::system_clock::now();
	//time = std::chrono::system_clock::to_time_t(now);
	//nowTime = std::localtime(&time);
	//errorWriter.open("./markerMatchError.txt", std::ios::app);
	//errorWriter << "[" << std::setw(2) << nowTime->tm_hour << ":"
	//	<< std::setw(2) << nowTime->tm_min << ":"
	//	<< std::setw(2) << nowTime->tm_sec << "]\n";
	//errorWriter << "number\t" << "source pt\t\t\t" << "target pt\t\t\t\t" << "error" << std::endl;

	//for (auto i = 0; i < sourceMatchedPoints.size(); ++i)
	//{
	//	auto error = glm::distance(glm::dvec3(transformation *
	//		glm::dvec4(sourceMatchedPoints[i], 1.0)), targetMatchedPoints[i]);
	//	printf("%d: %f\n", i, error);

	//	//2023/06/08 WYX
	//	errorWriter << i << "\t(" << (transformation * glm::dvec4(sourceMatchedPoints[i], 1.0))[0] << ", "
	//		<< (transformation * glm::dvec4(sourceMatchedPoints[i], 1.0))[1] << ", "
	//		<< (transformation * glm::dvec4(sourceMatchedPoints[i], 1.0))[2] << ")\t("
	//		<< targetMatchedPoints[i].x << ", " << targetMatchedPoints[i].y << ", " << targetMatchedPoints[i].z << ")\t\t"
	//		<< error << std::endl;
	//}
	//errorWriter << std::endl;
	//errorWriter.close();
	////==================== debug over ====================

	//应该去掉最大误差再配准一次



	return transformation;
}

//求出转换矩阵
glm::dmat4 estimateRigidTransformation(const int& matchedNum,
	const std::vector<glm::dvec3>& sourcePoints, const std::vector<glm::dvec3>& targetPoints)
{
	Eigen::Matrix<double, 4, 4> C1 = Eigen::Matrix<double, 4, 4>::Zero();
	Eigen::Matrix<double, 4, 4> C2 = Eigen::Matrix<double, 4, 4>::Zero();//初始化两个零矩阵
	auto c1 = C1.data();
	auto c2 = C2.data();

	glm::dvec3 a;
	glm::dvec3 b;
	for (auto i = 0; i < matchedNum; ++i)
	{
		a = sourcePoints[i];
		b = targetPoints[i];
		double axbx = a.x * b.x;
		double ayby = a.y * b.y;
		double azbz = a.z * b.z;
		double axby = a.x * b.y;
		double aybx = a.y * b.x;
		double axbz = a.x * b.z;
		double azbx = a.z * b.x;
		double aybz = a.y * b.z;
		double azby = a.z * b.y;
		c1[0] += axbx - azbz - ayby;
		c1[5] += ayby - azbz - axbx;
		c1[10] += azbz - axbx - ayby;
		c1[15] += axbx + ayby + azbz;
		c1[1] += axby + aybx;
		c1[2] += axbz + azbx;
		c1[3] += aybz - azby;
		c1[6] += azby + aybz;
		c1[7] += azbx - axbz;
		c1[11] += axby - aybx;

		c2[1] += a.z + b.z;
		c2[2] -= a.y + b.y;
		c2[3] += a.x - b.x;
		c2[6] += a.x + b.x;
		c2[7] += a.y - b.y;
		c2[11] += a.z - b.z;
	}

	c1[4] = c1[1];
	c1[8] = c1[2];
	c1[9] = c1[6];
	c1[12] = c1[3];
	c1[13] = c1[7];
	c1[14] = c1[11];
	c2[4] = -c2[1];
	c2[8] = -c2[2];
	c2[12] = -c2[3];
	c2[9] = -c2[6];
	c2[13] = -c2[7];
	c2[14] = -c2[11];

	C1 *= -2.0;
	C2 *= 2.0;

	const Eigen::Matrix<double, 4, 4> A = (0.25f / double(matchedNum)) * C2.transpose() * C2 - C1;
	
	const Eigen::EigenSolver<Eigen::Matrix<double, 4, 4>> es(A);
	//使用A矩阵	
	ptrdiff_t i;
	es.eigenvalues().real().maxCoeff(&i);
	const Eigen::Matrix<double, 4, 1> qmat = es.eigenvectors().col(i).real();
	const Eigen::Matrix<double, 4, 1> smat = -(0.5f / double(matchedNum)) * C2 * qmat;

	//构造元组 q s
	const Eigen::Quaternion<double> q(qmat(3), qmat(0), qmat(1), qmat(2));
	const Eigen::Quaternion<double> s(smat(3), smat(0), smat(1), smat(2));

	const Eigen::Quaternion<double> t = s * q.conjugate();

	const Eigen::Matrix<double, 3, 3> R(q.toRotationMatrix());

	//得出转换矩阵
	glm::dmat4 transformation;
		for (auto i = 0; i < 3; ++i)
		for (auto j = 0; j < 3; ++j)
			transformation[j][i] = R(i, j);

	transformation[3][0] = -t.x();
	transformation[3][1] = -t.y();
	transformation[3][2] = -t.z();

	return transformation;
}

//检测图像中的标记物，将检测到的坐标和半径保存在markers向量中
void detectMarkers(std::vector<std::pair<glm::dvec2, float>>& markers, 	const cv::Mat& image)
{
	//需要cv::Mat图像作为输入
	markers.clear();            //清除 markers矩阵，是存储了坐标和半径信息的向量
	findCircles(markers, image);//查找圆形标记物，将信息存到markers中
}

//检测图像中的圆圈
void findCircles(std::vector<std::pair<glm::dvec2, float>>& markers, const cv::Mat& image)
{
	//auto gamma_Transform = [](const cv::Mat& inputMat, cv::Mat& ouputMat, const float& factor)
	//{
	//	unsigned char lUT[256];
	//	for (int i = 0; i < 256; i++)
	//	{
	//		float f = (i + 0.5f) / 255.0f;
	//		f = (float)(pow(f, factor));
	//		lUT[i] = cv::saturate_cast<uchar>(f * 255.0f - 0.5f);
	//	}

	//	cv::Mat tem = inputMat.clone();
	//	cv::MatIterator_<uchar> iterator = tem.begin<uchar>();
	//	cv::MatIterator_<uchar> iteratorEnd = tem.end<uchar>();
	//	for (; iterator != iteratorEnd; iterator++)
	//	{
	//		*iterator = lUT[*iterator];
	//	}

	//	ouputMat = tem.clone();
	//};
	auto width = image.cols;
	auto height = image.rows;

	//高斯滤波
	cv::Mat smoothedImage = image.clone();
	//gamma_Transform(smoothedImage, smoothedImage, 0.8f);
	if (1)
	{
		cv::GaussianBlur(smoothedImage, smoothedImage, cv::Size(3, 3), 0.0);
	}

	//二值化图像
	cv::Mat mask;
	threshold(smoothedImage, mask, 50, 255, cv::THRESH_BINARY);

	//开运算
	//cv::morphologyEx(mask, mask, cv::MORPH_OPEN, cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3, 3)));//8
	//cv::adaptiveThreshold(smoothedImage, mask, 255, cv::AdaptiveThresholdTypes::ADAPTIVE_THRESH_GAUSSIAN_C, cv::THRESH_BINARY, 5, 0.0);

	cv::Mat edge; 	//原始图像边缘
	cv::Mat maskEdge; //二值化图像边缘
	Canny(smoothedImage, edge, 30, 90, 3);
	Canny(mask, maskEdge, 30, 90, 3);

	// 高斯变换
	double guassdeta = 1;		// 0.6
	cv::Mat k[5];
	for (int i = 0; i < 5; i++)
		k[i] = cv::Mat::zeros(height, width, CV_32FC1);

	cv::Mat imgsobel = cv::Mat::zeros(height, width, CV_32FC1);
	Gradient(smoothedImage, imgsobel);

	convolve_gauss(imgsobel, k[0], guassdeta, DERIV_R);
	convolve_gauss(imgsobel, k[1], guassdeta, DERIV_C);
	convolve_gauss(imgsobel, k[2], guassdeta, DERIV_RR);
	convolve_gauss(imgsobel, k[3], guassdeta, DERIV_RC);
	convolve_gauss(imgsobel, k[4], guassdeta, DERIV_CC);


	std::vector<std::vector<cv::Point>> contours; //原始图像轮廓
	std::vector<std::vector<cv::Point>> maskContours; //二值化图像轮廓
	//cv::findContours(edge, contours, cv::RETR_LIST, cv::CHAIN_APPROX_NONE);
	std::vector<cv::Vec4i> hierachy;
	//cv::findContours(edge, contours, hierachy, cv::RETR_TREE, cv::CHAIN_APPROX_NONE);
	//wxm对照calib
	cv::findContours(edge, contours, hierachy, cv::RETR_CCOMP, cv::CHAIN_APPROX_NONE);
	std::vector<std::vector<cv::Point>> contoursFiltered;

	//筛选没有子轮廓的轮廓
	for (int i = 0; i < contours.size(); i++)
	{
		if (hierachy[i][2] == -1)
			contoursFiltered.push_back(contours[i]);
	}

	cv::findContours(maskEdge, maskContours, cv::RETR_LIST, cv::CHAIN_APPROX_NONE);

	//xy方向的一阶导
	//cv::Mat sobelX;
	//cv::Mat sobelY;
	//cv::Sobel(image, sobelX, CV_32FC1, 1, 0, 1);
	//cv::Sobel(image, sobelY, CV_32FC1, 0, 1, 1);

	std::vector<cv::RotatedRect> circles;
	//亚像素边缘
	std::vector<std::vector<cv::Point2f>> subpixelContours;//mzy 20230923 临时改回Point，之后要改成Point2f
	std::vector<std::pair<glm::dvec2, float>> markersTmp;

	for (const auto& contour : contoursFiltered)
	{
		//for (cv::Point p : contour)
		//{
		//	if (p.x > 2480 && p.x < 2502 && p.y > 1647 && p.y < 1678)
		//		bool debugBegin = true;
		//}

		if (contour.size() >= 5)
		{
			//auto rect = cv::boundingRect(contour);
			auto rRect = cv::minAreaRect(contour).size; //最小外接矩形
			auto rectTh = 8;//13->4

			if (rRect.width > rectTh || rRect.height > rectTh)//WYX 2023/05/25 8 20(little) 10(baqiuMarker, strict)
			{
				//if (rRect.width < rectTh || rRect.height < rectTh)
				//{
				//	continue;
				//}

				//逐个轮廓
				auto area = fabs(cv::contourArea(contour)); //面积
				auto length = cv::arcLength(contour, true); //周长
				float diameter = 0; //直径
				for (const cv::Point& p1 : contour)
				{
					for (const cv::Point& p2 : contour)
					{
						float dis = norm(p1 - p2);
						if (dis > diameter)
						{
							diameter = dis;
						}
					}
				}
				//应该用拟合误差
				//if ( (area / length > length / SLS::pi / 4.0f * 0.87f && isCaptureBaqiuMarkers) 
				//	|| (area / length > length / SLS::pi / 4.0f * 0.885f && !isCaptureBaqiuMarkers)) //2023/06/07 WYX 圆形度约束 原来是0.6f  2023/08/16 0.87f
				

				//圆度约束0.65
				if (area / (M_PI * diameter * diameter / 4) > 0.65)
				{
					auto circle = cv::fitEllipse(contour);
					glm::ivec2 center(circle.center.x + 0.5f, circle.center.y + 0.5f);
					if (center.x >= 0 && center.y >= 0 && center.x < width && center.y < height)
					{
#pragma region wxm1
						//// 判据: 计算凸包 凸包的长度不应小于6
						//std::vector<cv::Point> hull;
						//cv::convexHull(contour, hull, false);
						//if (hull.size() < 6)
						//{
						//	continue;
						//}

						//std::vector<std::vector<cv::Point>> debugVec(1, hull);
						//cv::Mat debug1(2160, 4096, CV_8UC1);
						//drawContours(debug1, debugVec, -1, cv::Scalar::all(255), 0);
						//cv::RotatedRect minEllipse_hull = cv::fitEllipse(cv::Mat(hull));
						//double hull_length = cv::arcLength(hull, true);
						//double hull_area = cv::contourArea(hull);

						//auto distance = [](double x1, double y1, double x2, double y2) -> double
						//{
						//	double dx = x2 - x1;
						//	double dy = y2 - y1;
						//	return (sqrt(dx * dx + dy * dy));
						//};

						//auto fitError = [&distance](const std::vector<cv::Point>& pts, const cv::RotatedRect& ellipse) -> double
						//{
						//	const double& centerX = ellipse.center.x;
						//	const double& centerY = ellipse.center.y;

						//	// 椭圆半焦距
						//	const double& c = sqrt(fabs(ellipse.size.width * ellipse.size.width
						//		- ellipse.size.height * ellipse.size.height)) / 2.0;
						//	double angle = ellipse.angle;
						//	if (angle > 180) {
						//		angle = angle - 180;
						//	}
						//	if (ellipse.size.width < ellipse.size.height)
						//		angle = angle - 90;
						//	angle = angle * CV_PI / 180;

						//	// 计算两个焦点的位置
						//	const double& xc1 = centerX - c * cos(angle);
						//	const double& yc1 = centerY - c * sin(angle);
						//	const double& xc2 = centerX + c * cos(angle);
						//	const double& yc2 = centerY + c * sin(angle);
						//	double totalError = 0;

						//	std::vector<double> err_vec;						// 算标准差的时候 那些误差比较好的点就不要往里面放了
						//	for (int i = 0; i < (int)pts.size(); i++) {
						//		double err = fabs(distance(xc1, yc1, pts[i].x, pts[i].y)
						//			+ distance(xc2, yc2, pts[i].x, pts[i].y) - ellipse.size.height);
						//		totalError += pow(err, 2.0);
						//		if (err > 0.17)
						//			err_vec.push_back(err);
						//	}

						//	return (sqrt(totalError / double(pts.size())));
						//};

						//// 判据：椭圆拟合误差过大的不要
						//double aveError = fitError(contour, circle);
						//if (aveError > 1.0f)
						//{
						//	continue;
						//}

						//// 判据：拟合椭圆圆心过分靠近图像边缘的不要
						//double sideLength = circle.size.height * 1.1;
						////if (measuredTargetState != noCoded){
						//if (circle.center.x - 0.5 * sideLength < 0
						//	|| circle.center.y - 0.5 * sideLength < 0
						//	|| circle.center.x + 0.5 * sideLength > width
						//	|| circle.center.y + 0.5 * sideLength > height)
						//{
						//	continue;
						//}

						//// 判据：凸包拟合椭圆的面积和凸包的面积之差不应超过20%
						//if (1)
						//{
						//	if (fabs(hull_area - CV_PI * minEllipse_hull.size.width * minEllipse_hull.size.height / 4.0) / hull_area > 0.3)
						//	{
						//		continue;
						//	}
						//}


						//auto calFarthestDist2Hull = [](std::vector<cv::Point> contour, std::vector<cv::Point>& hull) -> double
						//{
						//	if (hull.empty())
						//		cv::convexHull(contour, hull, false);

						//	std::vector<double> dist(contour.size());
						//	double maxDist = 0;

						//	std::vector<cv::Point>::iterator it_c = contour.begin();
						//	for (int i = 0; i < contour.size(); i++)
						//	{
						//		dist[i] = cv::pointPolygonTest(hull, cv::Point2f(*it_c++), true);
						//		if (dist[i] > maxDist)
						//		{
						//			maxDist = dist[i];
						//		}
						//	}

						//	return maxDist;
						//};
						//// 判据:凸包拟合后 轮廓上的点距凸包的最大距离不应超过拟合椭圆短轴的1/15
						//if (1)
						//{
						//	if (calFarthestDist2Hull(contour, hull) > circle.size.height / 10.0)
						//		//if (calFarthestDist2Hull(contour, hull) > 2) //2023/08/17 JZH
						//	{
						//		continue;
						//	}
						//}
#pragma endregion wxm2
						//中心约束
						if (mask.at<uchar>(center.y, center.x) >= 0)//==255
						{
							//计算contour内部灰度均值，小的不要
							double meanGrayTh = 140; //2023/08/17 180
							cv::Mat contourMask = cv::Mat::zeros(mask.size(), mask.type());
							std::vector<std::vector<cv::Point>> contoursTmp;
							contoursTmp.push_back(contour);
							cv::drawContours(contourMask, contoursTmp, 0, cv::Scalar::all(255), -1, 8);
							cv::Scalar meanGray = cv::mean(mask, contourMask);
							double innerMeanGray = meanGray[0];
							OutputDebugStringA(("innerMeanGray = " + std::to_string(innerMeanGray) + "\r\n").c_str());
							if (innerMeanGray < meanGrayTh)
							{
								continue;
							}

							auto inside = true;
							std::vector<cv::Point2f> subpixelPoints;

							//亚像素边缘提取
							cv::Point2d centerPt;
							bool test = FindSubEdge(contour, centerPt, k, 1.0, 20, 2.0);
							
							for (const auto& point : contour)
							{
								subpixelPoints.push_back(cv::Point2f(point.x, point.y));
							}
							if (inside && subpixelPoints.size() >= 5)
							{
								subpixelContours.emplace_back();
								for (const auto& subpixelPoint : subpixelPoints)
									subpixelContours.back().push_back(subpixelPoint);
								//printf("(%f %f) ", circle.center.x, circle.center.y);
								circle = cv::fitEllipse(subpixelPoints);
								//printf("(%f %f)\n", circle.center.x, circle.center.y);
								circles.push_back(circle);
							}
							markersTmp.emplace_back(glm::dvec2(centerPt.x, centerPt.y), circle.size.area());

							////逐个点
							//for (const auto& point : contour)
							//{
							//	if (point.x < 2 || point.y < 2 ||
							//		point.x >= width - 2 || point.y >= height - 2)
							//	{
							//		inside = false;
							//		break;
							//	}
							//	//内部点
							//	else
							//	{
							//		//auto delta = subpixelFit(point, image,
							//		//	sobelX.at<float>(point.y, point.x), sobelY.at<float>(point.y, point.x));
							//		//亚像素拟合
							//		auto delta = subpixelFit(point, image, sobelX, sobelY);
							//		if (abs(delta.x) < 0.5 && abs(delta.y) < 0.5)
							//			subpixelPoints.push_back(cv::Point2f(point.x, point.y) + delta);
							//	}
							//}
							//if (inside && subpixelPoints.size() >= 5)
							//{
							//	//亚像素边缘
							//	subpixelContours.emplace_back();
							//	for (const auto& subpixelPoint : subpixelPoints)
							//		subpixelContours.back().push_back(subpixelPoint);
							//	//printf("(%f %f) ", circle.center.x, circle.center.y);
							//	circle = cv::fitEllipse(subpixelPoints);
							//	//printf("(%f %f)\n", circle.center.x, circle.center.y);
							//	circles.push_back(circle);
							//}
						}
					}
				}

			}
		}
	}


	std::vector<std::pair<cv::RotatedRect, float>> cleanCircles;
	cleanCircles.clear();
	std::vector<float> circleArea;
	circleArea.clear();
	std::vector<bool> processed(circles.size(), false);
	float _minArea = 1000.0;
	for (auto i = 0; i < circles.size(); ++i)
	{
		if (!processed[i])   //检测是否有圆未处理
		{
			auto minArea = circles[i].size.area();
			auto minAreaCircle = circles[i];
			std::pair<glm::dvec2, float> minAreaMarker = markersTmp[i];//初始化为当前圆
			for (auto j = 0; j < circles.size(); ++j)
			{
				if (i != j && !processed[j] && glm::distance(glm::vec2(circles[i].center.x,
					circles[i].center.y), glm::vec2(circles[j].center.x, circles[j].center.y)) < 5.0f)
				{
					if (circles[j].size.area() < minArea)
					{
						minArea = circles[j].size.area();
						minAreaCircle = circles[j];
						minAreaMarker = markersTmp[j];
					}
					processed[j] = true;
				}
			}
			if (minArea < _minArea)
			{
				_minArea = minArea;
			}
			cleanCircles.push_back(make_pair(minAreaCircle, minArea));
			markers.push_back(minAreaMarker);
			circleArea.emplace_back(minArea);
			processed[i] = true;
		}
	}



	if (1)          //debug
	{
		cv::Mat contourImage;    //绘制轮廓
		cv::cvtColor(image, contourImage, cv::COLOR_GRAY2RGB);

		cv::Mat circleImage;
		cv::cvtColor(image, circleImage, cv::COLOR_GRAY2RGB);

		for (auto i = 0; i < subpixelContours.size(); ++i)
		{
			std::vector<cv::Point2i> contourInt;
			for (auto& p : subpixelContours[i])
			{
				//检测坐标是否超出范围
				if (p.x < 0 || p.x>4095 || p.y < 0 || p.y>2159)
				{
					OutputDebugStringA(("p:[" + std::to_string(p.x) + ", " + std::to_string(p.y) + "]\r\n").c_str());
					break;
				}
				contourInt.push_back(static_cast<cv::Point2i>(p));
				//转换类型，添加到contourint向量中
				//OutputDebugStringA(("p:[" + std::to_string(p.x) + ", " + std::to_string(p.y) + "]\r\n").c_str());
			}
			std::vector<std::vector<cv::Point2i>> contourIntTmp;
			contourIntTmp.push_back(contourInt);
			cv::drawContours(contourImage, contourIntTmp, -1, cv::Scalar(rand() % 255,
				rand() % 255, rand() % 255), 1);
			//调用drawcontours函数
		}


		for (const auto& circle : cleanCircles)
		{
			//使用sllipse函数绘制椭圆
			cv::ellipse(circleImage, circle.first, cv::Scalar(0, 0, 255));
			//绘制小圆点
			cv::circle(circleImage, cv::Point(circle.first.center.x, circle.first.center.y), 1, cv::Scalar(255, 0, 0), -1);
		}
		imwrite(std::string("image.png"), smoothedImage);
		imwrite(std::string("mask.png"), mask);
		imwrite(std::string("contour.png"), contourImage);
		imwrite(std::string("edge.png"), edge);
		imwrite(std::string("circle.png"), circleImage);//保存
		cv::Mat markerImage(image.rows, image.cols, CV_8UC3, cv::Scalar(0, 0, 0));
		cv::cvtColor(smoothedImage, markerImage, cv::COLOR_GRAY2BGR);
		for (const auto& marker : markers)
			cv::circle(markerImage, cv::Point(marker.first.x, marker.first.y), 1, cv::Scalar(0, 0, 255), -1);//红色圆心
		imwrite(std::string("markerImg.png"),markerImage);
		}
}

//4.0之前用的亚像素边缘提取  （索贝尔算子）      返回cv::Point2f 亚像素二维坐标
cv::Point2f subpixelFit(const cv::Point& point, const cv::Mat& image,	const cv::Mat& sobelX, const cv::Mat& sobelY)
{
	auto sumWeightsX = 0.0;
	auto sumWeightsY = 0.0;
	auto sumX = 0.0;
	auto sumY = 0.0;
	auto differenceX = sobelX.at<float>(point.y, point.x);
	auto differenceY = sobelY.at<float>(point.y, point.x);   //获取像素差异值
	auto computeSum = [&sumX, &sumY, &sumWeightsX, &sumWeightsY, &sobelX, &sobelY, &point] //计算加权和
	(const int& x, const int& y)
	{
		sumX += sobelX.at<float>(y, x) * (x - point.x);
		sumY += sobelY.at<float>(y, x) * (y - point.y);
		sumWeightsX += sobelX.at<float>(y, x);
		sumWeightsY += sobelY.at<float>(y, x);
	};
	auto angle = glm::degrees(atan2f(differenceY, differenceX));  //转换为度数制
	if ((angle > 45.0f && angle < 135.0f) || (angle > -135.0f && angle < -45.0f))
		for (auto y = point.y - 2; y <= point.y + 2; ++y)
			computeSum(point.x, y);
	else
		for (auto x = point.x - 2; x <= point.x + 2; ++x)
			computeSum(x, point.y);
	return cv::Point2f(sumX / sumWeightsX, sumY / sumWeightsY);//返回通过亚像素拟合的坐标
}

//归一化 
glm::dvec2 normalize(const glm::dvec2&point, SLSLens lens)  //相机中心坐标cc，焦距坐标fc，畸变系数kc
{
	glm::dvec2 normalizedPoint;  //存储归一化坐标
	glm::dvec2 x_kk = point;   //输入坐标赋值给x_kk
	glm::dvec2 cc = lens.cc;
	glm::dvec2 fc = lens.fc;
	glm::dvec4 kc = lens.kc;
	normalizedPoint = undistort((x_kk - cc) / fc, kc);
	return normalizedPoint;
}

//去畸变√
glm::dvec2 undistort(const glm::dvec2& x_kk, const glm::dvec4& kc)//需校正的输入坐标，畸变系数
{
	if (glm::length(kc) == 0.0)
		return x_kk;
	auto pt = x_kk;
	glm::vec2 delta;  //存储校正量
	for (auto i = 0; i < 20; ++i)
	{
		auto xx = pt.x * pt.x;
		auto yy = pt.y * pt.y;
		auto r2 = xx + yy;
		delta.x = 2.0f * kc.z * pt.x * pt.y + kc.w * (r2 + 2.0f * xx);
		delta.y = kc.z * (r2 + 2.0f * yy) + 2.0f * kc.w * pt.x * pt.y;
		auto dr = 1.0f + kc.x * r2 + kc.y * r2 * r2;
		pt.x = (x_kk.x - delta.x) / dr;
		pt.y = (x_kk.y - delta.y) / dr;
		if (glm::distance(pt, x_kk) < 1e-8f)  //更新点坐标和输入坐标之间的距离是否小于阈值
			break;
	}
	return pt;
}

//极线匹配误差    K_inv表示逆相机内参矩阵向量
glm::dvec2 epipolarErr(const glm::dvec2& normalizedSourcePoint, 
	const glm::dvec2& normalizedMatchPoint, glm::dmat3 essentialMatrix,vector<glm::dmat3> K_inv) //计算本质矩阵
{
	glm::dmat3 fundamentalMatrix = glm::transpose(K_inv[1]) * essentialMatrix * K_inv[0];//计算两幅二维图像位置关系
	glm::dvec3 epipolarLine = fundamentalMatrix * glm::dvec3(normalizedSourcePoint, 1.0);
	//归一化匹配点坐标与其转置相乘，得到极线
	//cout << "\nfundamentalMatrix: " << endl;
	//cout << fundamentalMatrix[0][0] << "\t" << fundamentalMatrix[0][1] << "\t" << fundamentalMatrix[0][2] << endl;
	//cout << fundamentalMatrix[1][0] << "\t" << fundamentalMatrix[1][1] << "\t" << fundamentalMatrix[1][2] << endl;
	//cout << fundamentalMatrix[2][0] << "\t" << fundamentalMatrix[2][1] << "\t" << fundamentalMatrix[2][2] << endl;
	//cout << "\nepipolarLine: " << endl;
	//cout << epipolarLine.x << " " << epipolarLine.y << " " << epipolarLine.z << endl << endl;
	
	return glm::dvec2(fabs(glm::dot(glm::dvec3(normalizedMatchPoint, 1.0), epipolarLine)), sqrtf(epipolarLine.x * epipolarLine.x + epipolarLine.y * epipolarLine.y));
	//返回极线的误差分量
}

//计算极线误差
glm::dvec2 epipolarErr2(const glm::dvec2& normalizedSourcePoint, 
	const glm::dvec2& normalizedMatchPoint, glm::dvec3& epipolarLine, glm::dmat3 essentialMatrix, vector<glm::dmat3> K_inv)
{
	glm::dmat3 fundamentalMatrix = glm::transpose(K_inv[1]) * essentialMatrix * K_inv[0];//存储计算得到的基础矩阵
	epipolarLine = glm::transpose(fundamentalMatrix) * glm::dvec3(normalizedMatchPoint, 1.0);//计算极线
	return glm::dvec2(fabs(glm::dot(glm::dvec3(normalizedSourcePoint, 1.0), epipolarLine)), sqrtf(epipolarLine.x * epipolarLine.x + epipolarLine.y * epipolarLine.y));
	//源点坐标与极线点积的绝对值是误差的第一个分量，极线在xy防线上的平方和和平方根是第二分量
}

//检查极线约束    参数要有归一化源点坐标、匹配点坐标，极线约束的阈值，本质矩阵，逆相机内参矩阵的向量
bool checkEpipolar(const glm::dvec2& normalizedSourcePoint,
	const glm::dvec2& normalizedMatchPoint, const float& threshold, glm::dmat3 essentialMatrix, vector<glm::dmat3> K_inv)
{
	glm::dmat3 fundamentalMatrix = glm::transpose(K_inv[1]) * essentialMatrix * K_inv[0];
	glm::dvec3 epipolarLine = fundamentalMatrix * glm::dvec3(normalizedSourcePoint, 1.0);
	return fabs(glm::dot(glm::dvec3(normalizedMatchPoint, 1.0), epipolarLine)) /
		sqrt(epipolarLine[0] * epipolarLine[0] + epipolarLine[1] * epipolarLine[1]) < threshold;
	//点积绝对值除以极线的长度 
}

//检查标记点匹配
double checkMarkers(cv::Mat imageL, cv::Mat imageR, const glm::dvec2& sourceMarker,
	const vector<glm::dvec2>& targetMarkers, vector<double> targetEpipLs, vector<double> targetEpipRs, 
	glm::dvec2& rightMarker, double& matchEpipL, double& matchEpipR, const int& roiPara)
{
	//int roiPara = 40; //截取sourceMarkers为中心的区域，向外延伸roiPara个像素

	//判断roi是否在图像外
	if (sourceMarker.x - roiPara < 0 || sourceMarker.x + roiPara > imageL.cols
		|| sourceMarker.y - roiPara < 0 || sourceMarker.y + roiPara > imageL.rows)
	{
		rightMarker = glm::dvec2();
		return 0;
	}

	cv::Mat sourceROI = imageL(cv::Rect(sourceMarker.x - roiPara, sourceMarker.y - roiPara, 2 * roiPara, 2 * roiPara));
	cv::Scalar sourcegrayMean = cv::mean(sourceROI);//计算灰度平均值 
	double maxS = 0;
	for (auto i = 0; i < targetMarkers.size(); i++)
	{
		if (targetMarkers[i].x - roiPara < 0 || targetMarkers[i].x + roiPara > imageR.cols
			|| targetMarkers[i].y - roiPara < 0 || targetMarkers[i].y + roiPara > imageR.rows)
		{
			rightMarker = glm::dvec2();
			return 0;
		}

		cv::Mat targetROI = imageR(cv::Rect(targetMarkers[i].x - roiPara, targetMarkers[i].y - roiPara, 2 * roiPara, 2 * roiPara));
		cv::Scalar targetgrayMean = cv::mean(targetROI);
		double srcTmp = 0, tgtTmp = 0, numerator = 0, denominatorTmp1 = 0, denominatorTmp2 = 0;//numerator分子，denominator分母
		for (auto m = 0; m < targetROI.rows; m++)  
		{
			for (auto n = 0; n < targetROI.cols; n++)
			{
				srcTmp = sourceROI.at<uchar>(m, n) - sourcegrayMean[0];
				tgtTmp = targetROI.at<uchar>(m, n) - targetgrayMean[0];
				numerator += srcTmp * tgtTmp;
				denominatorTmp1 += pow(srcTmp, 2);
				denominatorTmp2 += pow(tgtTmp, 2);
				//计算源标记点区域和目标点的灰度差，进行相似度计算
			}
		}
		//像素   

		double S = numerator / sqrt(denominatorTmp1 * denominatorTmp2);
		if (S > maxS)
		{
			rightMarker = targetMarkers[i];
			matchEpipL = targetEpipLs[i];
			matchEpipR = targetEpipRs[i];
			maxS = S;
		}
	}
	//返回最大相似度
	return maxS;
}


//三角化源点和匹配点的空间坐标
glm::dvec3 triangulate(const glm::dvec2& sourcePoint,
	const glm::dvec2& matchedPoint, vector<SLSLens> lens)    //参数：源点坐标、匹配点坐标、镜头参数
{
	//auto kc = lens[0].kc;
	//if (!considerDistortion)
	//	kc = glm::dvec4();    镜头修正，如果不考虑畸变，kc设置为0
	auto xt = glm::dvec3(normalize(glm::dvec2(sourcePoint), lens[0]), 1.0);
	//kc = lens[1].kc;
	//if (!considerDistortion)
	//	kc = glm::dvec4();
	auto xtt = glm::dvec3(normalize(glm::dvec2(matchedPoint), lens[1]), 1.0);//转换为相机坐标系下的坐标
	
	auto u = lens[1].R * xt;
	auto n_xt2 = glm::dot(xt, xt);
	auto n_xtt2 = glm::dot(xtt, xtt);

	auto DD = n_xt2 * n_xtt2 - glm::pow(glm::dot(u, xtt), 2.0);
	auto dot_uT = glm::dot(u, lens[1].t);
	auto dot_xttT = glm::dot(xtt, lens[1].t);
	auto dot_xttu = glm::dot(u, xtt);    //向量点积
	auto NN1 = dot_xttu * dot_xttT - n_xtt2 * dot_uT;
	auto NN2 = n_xt2 * dot_xttT - dot_uT * dot_xttu;
	auto X1 = NN1 / DD * xt;
	auto X2 = glm::transpose(lens[1].R) * (NN2 / DD * xtt - lens[1].t);
	//cout << "DD: " << DD<< endl;
	//cout << "NN1: " << NN1 << endl;
	//cout << "NN2: " << NN2 << endl;
	//cout << "X1: " << X1.x << " " << X1.y << " " << X1.z << endl;
	//cout << "X2: " << X2.x << " " << X2.y << " " << X2.z << endl;

	return 0.5 * (X1 + X2);   //返回源点和匹配点的三角化结果的平均值
}













//拟合标定平面
bool fitMarkerPlane(const glm::dvec3& basePoint,  const vector<glm::dvec3>& points, 
	const int& maxIterations, const double& minInnerRatio, const double& maxError,
	glm::dvec3& planePoint, glm::dvec3& planeNormal)
{
	if (points.size() < 3)
	{

		planePoint = glm::dvec3();
		planeNormal = glm::dvec3();
		return false;
	}

	auto maxInnnersNum = 0;
	int randomNum[3];
	int tempRandomNum[3];
	glm::dvec3 innerPlaneNormal;
	glm::dvec3 innerPlanePoint;
	std::vector<glm::dvec3> selectedPoints(3);
	selectedPoints[0] = basePoint;
	for (int i = 0; i < points.size(); ++i)
	{
		for (int j = i + 1; j < points.size(); ++j)
		{
			if (points[i] == basePoint || points[j] == basePoint)
				continue;
			selectedPoints[1] = points[i];
			selectedPoints[2] = points[j];

			std::pair<glm::dvec3, glm::dvec3> planeCoeff_temp;
			//glm::dvec3 planePoint_temp, planeNormal_temp;
			planeCoeff_temp = fitPlane(selectedPoints);
			planePoint = planeCoeff_temp.first;
			planeNormal = planeCoeff_temp.second;

			auto innersNum = 0;

			for (int k = 0; k < points.size(); ++k)
			{
				auto distance = glm::dot(points[k] - planePoint, planeNormal);
				if (fabs(pointPlaneDistance(points[k], planePoint, planeNormal)) < maxError)
					++innersNum;
			}
			if (innersNum > maxInnnersNum)
			{
				maxInnnersNum = innersNum;
				innerPlaneNormal = planeNormal;
				innerPlanePoint = planePoint;
			}
		}
	}
	/*for (auto i = 0; i < maxIterations; ++i)
	{
		auto innersNum = 0;
		tempRandomNum[0] = -1;
		tempRandomNum[1] = -1;
		for (auto j = 0; j < 3; ++j)
		{
			randomNum[j] = rand() % points.size();
			while (randomNum[j] == tempRandomNum[0] || randomNum[j] == tempRandomNum[1])
				randomNum[j] = rand() % points.size();
			tempRandomNum[j] = randomNum[j];
		}

		for (auto j = 0; j < 3; ++j)
			selectedPoints[j] = points[randomNum[j]];
		selectedPoints[0] = basePoint;

		auto [planePoint, planeNormal] = fitPlane(selectedPoints);

		for (int j = 0; j < points.size(); ++j)
		{
			auto distance = glm::dot(points[j] - planePoint, planeNormal);
			if (fabs(pointPlaneDistance(points[j], planePoint, planeNormal)) < maxError)
				++innersNum;
		}
		if (innersNum > maxInnnersNum)
		{
			maxInnnersNum = innersNum;
			innerPlaneNormal = planeNormal;
			innerPlanePoint = planePoint;
		}
		if ((float)innersNum / points.size() > minInnerRatio)
			break;
	}*/

	//
	std::ofstream debug("./debug4PlaneFitting.txt", std::ios::ate);
	debug << "maxError: " << maxError << std::endl;
	debug << "planeNormal: " << innerPlaneNormal.x << " " << innerPlaneNormal.y << " " << innerPlaneNormal.z << std::endl;
	debug << "innerPlanePoint: " << innerPlanePoint.x << " " << innerPlanePoint.y << " " << innerPlanePoint.z << std::endl;

	std::vector<glm::dvec3> innerPoints;
	for (const auto& point : points)
	{
		if (fabs(pointPlaneDistance(point, innerPlanePoint, innerPlaneNormal)) < maxError)
		{
			innerPoints.push_back(point);
			debug << "point: " << point.x << " " << point.y << " " << point.z << std::endl;
			debug << "pointPlaneDistance: " << fabs(pointPlaneDistance(point, innerPlanePoint, innerPlaneNormal)) << std::endl;
		}

	}
	debug << "innerPoints.size=" << innerPoints.size() << std::endl;
	debug.close();


	printf("innersNum %zd\n", innerPoints.size());

	//大小如果<6，表示内点数量不足
	if (innerPoints.size() >= 6)//20240129 3->6 YDR
	{
		std::pair<glm::dvec3, glm::dvec3> planeCoeff;
		planeCoeff = fitPlane(innerPoints);
		planePoint = planeCoeff.first;
		planeNormal = planeCoeff.second;
		//平面拟合 得到平面系数
		return true;
	}
	else
	{
		planePoint = glm::dvec3();
		planeNormal = glm::dvec3();
		return false;
	}
}

std::pair<glm::dvec3, glm::dvec3> fitPlane(const std::vector<glm::dvec3>& points)
{
	auto parasNum = 3;
	if (points.size() >= parasNum)
	{
		glm::dvec3 normal;
		glm::dvec3 center;  //点云质心
		for (auto i = 0; i < points.size(); ++i)
			center = i / (i + 1.0) * center + 1.0 / (i + 1.0) *
			glm::dvec3(points[i].x, points[i].y, points[i].z);
		Eigen::MatrixXd matrix(points.size(), 3);
		for (auto i = 0; i < points.size(); ++i)
		{
			matrix(i, 0) = points[i].x - center.x;
			matrix(i, 1) = points[i].y - center.y;
			matrix(i, 2) = points[i].z - center.z;
		}

		// 对矩阵进行特征值分解并且提取特征向量
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(matrix.transpose() * matrix);
		auto eigenVectors = eigen.eigenvectors();
		normal.x = eigenVectors(0, 0);
		normal.y = eigenVectors(1, 0);
		normal.z = eigenVectors(2, 0);
		//cout << "======平面参数：" << center.x << " " << center.y << " " << center.z << " "
			//<< normal.x << " " << normal.y << " " << normal.z << endl;
		return std::make_pair(glm::dvec3(center.x, center.y, center.z),
			glm::normalize(normal));
		//返回法向量和质心
	}
	else
		return std::make_pair(glm::dvec3(), glm::dvec3());
}

//计算点到平面的距离
double pointPlaneDistance(const glm::dvec3& point,
	const glm::dvec3& planePoint, const glm::dvec3& planeNormal)
{
	return fabs(glm::dot(point - planePoint, planeNormal));
}
//点坐标-平面点坐标得到向量，点积之后得到投影长度

//矫正标记点的中心坐标
std::string correctMarkerCenters(std::vector<glm::dvec2>& markers,
	const std::vector<glm::dvec3>& markerPoints, const glm::dvec3& planeNormal,
	const float& markerRadius, const SLSLens& lens)
{
	if (markers.size() != markerPoints.size())
		return "标志点中心与标志点坐标数量不同！\n";
	//检查大小是否一样

	auto samplesNum = 100;
	std::vector<glm::dvec3> circlePoints(samplesNum);
	for (auto i = 0; i < circlePoints.size(); ++i)
	{
		circlePoints[i].x = markerRadius * cosf((2.0f * i / samplesNum) * M_PI);
		circlePoints[i].y = markerRadius * sinf((2.0f * i / samplesNum) * M_PI);
	}

	//计算旋转矩阵
	auto rotateMatrix = glm::dmat3(glm::rotate(glm::dmat4(),
		glm::degrees(acos(glm::dot(planeNormal, glm::dvec3(0.0, 0.0, 1.0)))),
		glm::cross(glm::dvec3(0.0, 0.0, 1.0), planeNormal)));

	std::vector<cv::Point2f> projectedCirclePoints(samplesNum);
	for (auto i = 0; i < markers.size(); ++i)
	{
		auto projectedCircleCenter = project(markerPoints[i], lens);//计算投影中心
		for (auto j = 0; j < samplesNum; ++j)
		{
			auto projectedCirclePoint = project(rotateMatrix * circlePoints[j] +
				markerPoints[i], lens);
			projectedCirclePoints[j].x = projectedCirclePoint.x;
			projectedCirclePoints[j].y = projectedCirclePoint.y;
		}

		auto circle = cv::fitEllipse(projectedCirclePoints);//拟合椭圆
		markers[i].x += projectedCircleCenter.x - circle.center.x;
		markers[i].y += projectedCircleCenter.y - circle.center.y;//校正
		//printf("(%f %f) (%f %f)\n", markers[i].x, markers[i].y,
		//	projectedCircleCenter.x - circle.center.x,
		//	projectedCircleCenter.y - circle.center.y);
	}

	return "";
}

//3D投影到2D平面上
glm::dvec2 project(const glm::dvec3& point, const SLSLens& lens)
{
	auto transformedPoint = lens.R * point + lens.t; //得到变换后的点坐标
	auto xx = glm::dvec2(transformedPoint.x / transformedPoint.z,
		transformedPoint.y / transformedPoint.z);   //得到归一化二维坐标
	auto r2 = glm::dot(xx, xx);
	auto xxx = glm::dvec2(xx.x * (1.0f + lens.kc.x * r2 + lens.kc.y * r2 * r2) +
		2.0f * lens.kc.z * xx.x * xx.y + lens.kc.w * (r2 + 2.0f * xx.x * xx.x),
		xx.y * (1.0f + lens.kc.x * r2 + lens.kc.y * r2 * r2) +
		2.0f * lens.kc.w * xx.x * xx.y + lens.kc.z * (r2 + 2.0f * xx.y * xx.y));
	return lens.fc * xxx + lens.cc; //归一化坐标*焦距+主点偏移量        得到最终投影坐标
}

