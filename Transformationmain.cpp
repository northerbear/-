#define _CRT_SECURE_NO_WARNINGS 1 

#include "4.0Function.h"
#include "calibFunction.h"
#include "fitFunction.h"  
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace cv;


int	main()
{

	//---------------------------------------------------------相机参数---------------------------------------------------------
	vector<SLSLens> lens;
	getLens(lens);

	//读图
	string filenameL = "globalMarkerImgL.png";
	string filenameR = "globalMarkerImgR.png";
	// filenameL = "L1001.BMP";
	//string filenameR = "R1001.BMP";

	cv::Mat imageL = cv::imread(filenameL, 0);
	cv::Mat imageR = cv::imread(filenameR, 0);


	//---------------------------------------------------------重建三维点---------------------------------------------------------
	vector<glm::dvec3> markerPoints;
	vector<glm::dvec2> matchMarkersL, matchMarkersR;
	vector<double> epipL, epipR;
	vector<glm::dvec4>Transformation;

	markerPoints = reconstruct3dPoints(lens, imageL, imageR,
		matchMarkersL, matchMarkersR, epipL, epipR);

	return 0;
}

int main()
{

}