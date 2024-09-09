#pragma once
//sobel 卷积求圆边缘轮廓
#include "calibFunction.h"

using namespace std;
using namespace cv;

//求梯度图像
bool Gradient(cv::Mat src, cv::Mat& dst)
{
	cv::Mat img;
	src.convertTo(img, CV_32FC1);
	int w = src.cols;
	int h = src.rows;
	dst = cv::Mat::zeros(src.size(), CV_32FC1);

#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 1; i < h - 1; i++)
	{
		float* s = img.ptr<float>(i);
		float* s_up = img.ptr<float>(i - 1);
		float* s_dw = img.ptr<float>(i + 1);
		float* pr = dst.ptr<float>(i);
		for (int j = 1; j < w - 1; j++)
		{
			float sx = -s_up[j - 1] + s_up[j + 1] - 2.0 * s[j - 1] + 2.0 * s[j + 1] - s_dw[j - 1] + s_dw[j + 1];
			float sy = s_up[j - 1] + 2.0 * s_up[j] + s_up[j + 1] - s_dw[j - 1] - 2.0 * s_dw[j] - s_dw[j + 1];
			pr[j] = sqrt(sx * sx + sy * sy);
		}
	}
	return true;
}

// 卷积 高斯核
void convolve_gauss(cv::Mat src, cv::Mat& dst, double sigma, long deriv_type) {
	cv::Mat orders1;
	switch (deriv_type) {
	case DERIV_R:
		orders1 = (cv::Mat_<double>(1, 2) << 0, 1); break;
	case DERIV_C:
		orders1 = (cv::Mat_<double>(1, 2) << 1, 0); break;
	case DERIV_RR:
		orders1 = (cv::Mat_<double>(1, 2) << 0, 2); break;
	case DERIV_RC:
		orders1 = (cv::Mat_<double>(1, 2) << 1, 1); break;
	case DERIV_CC:
		orders1 = (cv::Mat_<double>(1, 2) << 2, 0); break;
	}
	gfilter(src, dst, sigma, orders1);
}


// Gaussian filtering and Gaussian derivative filters
void gfilter(cv::Mat src, cv::Mat& dst, double sigma, cv::Mat orders)
{
	cv::Mat LL;
	int dims = 2;
	cv::Mat kk;
	int sze = 0;

	//先是列，后是行

	for (int i = 0; i < dims; i++)
	{
		// 返回k为列向量
		sze = gaussiankernel((float)sigma, (int)orders.ptr<double>(0)[i], kk);
		//std::cout << "sze = " << sze << std::endl;

		// shift the dimension of the kernel
		if (i == 0)
		{
			//std::cout << "kk = " << kk << std::endl;
			filter2D(src, LL, src.depth(), kk, cv::Point(-1, -1), 0, 1);
		}
		else if (i == 1)
		{
			cv::Mat kk1(1, 2 * sze + 1, CV_64FC1);
			memcpy((double*)kk1.data, (double*)kk.data, sizeof(double) * (2 * sze + 1));
			//std::cout << "kk1 = " << kk1 << std::endl;
			filter2D(LL, dst, CV_32FC1, kk1, cv::Point(-1, -1), 0, 1);
		}
	}
}


// GAUSSIAN DERIVATIVE KERNEL - creates a gaussian deivative kernel.
int gaussiankernel(float sigma, int order, cv::Mat& outputArray)
{
	float sigma2 = sigma * sigma;
	float sigma2szeRatio = (float)(3.0f + 0.25 * order - 2.5 / (pow((double)(order - 6), 2.0) + pow((double)(order - 9), 2.0)));

	// calculate kernel size
	int	sz = int(ceil(float(sigma2szeRatio * sigma)));
	outputArray = cv::Mat::zeros(2 * sz + 1, 1, CV_64FC1);

	double* temp = new double[2 * sz + 1];
	double* temp1 = new double[2 * sz + 1];
	double* temp2 = new double[2 * sz + 1];
	double* tempk = outputArray.ptr<double>(0);
	double* temppart = new double[2 * sz + 1];
	for (int i = 0; i < (2 * sz + 1); i++) {
		temp[i] = -sz + i;	//	t[i] = -sze + i, t的范围从-sze到sze
	}

	// CALCULATE GAUSSIAN
	for (int i = 0; i < (2 * sz + 1); i++) {
		temp1[i] = exp(-(temp[i]) * (temp[i]) / (2.0f * sigma2));
		temp2[i] = temp[i] / (sigma * sqrt(2.0f));
	}

	switch (order) {
	case 0:
		for (int i = 0; i < (2 * sz + 1); i++) temppart[i] = 1.0; break;
	case 1:
		for (int i = 0; i < (2 * sz + 1); i++)
			temppart[i] = temp2[i] * 2;
		break;
	case 2:
		for (int i = 0; i < (2 * sz + 1); i++)
			temppart[i] = temp2[i] * temp2[i] * 4.0 - 2.0;
		break;
	default:
	{
		std::cout << "There is a problem!" << std::endl;
		break;
	}
	}

	// apply Hermite polynomial to gauss
	for (int i = 0; i < (2 * sz + 1); i++)
		tempk[i] = (pow(-1.0, (double)order)) * (temppart[i] * temp1[i]);
	// Normalize
	double Sum = 0;
	for (int i = 0; i < (2 * sz + 1); i++)
	{
		Sum += temp1[i];
	}

	double norm_default = 1.0 / Sum;
	double norm_hermite = 1.0 / (pow((sigma * sqrt(2.0f)), order));
	for (int i = 0; i < (2 * sz + 1); i++)
	{
		tempk[i] = tempk[i] * (norm_default * norm_hermite);
	}

	delete[] temp;
	delete[] temp1;
	delete[] temp2;
	tempk = NULL;
	delete[] temppart;

	return sz;
}

bool FindSubEdge(vector<cv::Point> pts, cv::Point2d& centerPt, cv::Mat* k,
	double FitErrorThreshold, double m_minNumpt, double T)
{
	int count = (int)pts.size();

	vector<cv::Point2f> d2;
	vector<cv::Point2f> dTmp(count);
	vector<int> idTmp(count);

#pragma omp parallel for
	for (int m = 0; m < count; m++)
	{
		idTmp[m] = -1;
		int j = pts[m].x;
		int i = pts[m].y;
		//printf("pos: %d %d\n", j, i);

		int offset = i * k[0].cols + j;

		double rx = k[0].ptr<float>(0)[offset];
		double ry = k[1].ptr<float>(0)[offset];
		double rxx = k[2].ptr<float>(0)[offset];
		double rxy = k[3].ptr<float>(0)[offset];
		double ryy = k[4].ptr<float>(0)[offset];

		cv::Mat rrr = (cv::Mat_<double>(2, 2) << rxx, rxy, rxy, ryy);
		cv::Mat eigs, eigVecs;
		eigen(rrr, eigs, eigVecs);
		//cout << "rrr\n" << rrr << endl;
		//cout << "tezhengzhi\n" << eigs << endl;

		double nx = 0, ny = 0, eigValue = 0;
		if (fabs(eigs.ptr<double>(0)[0]) > fabs(eigs.ptr<double>(0)[1]))
		{
			nx = eigVecs.ptr<double>(0)[0];
			ny = eigVecs.ptr<double>(0)[1];
			eigValue = eigs.ptr<double>(0)[0];
		}
		else
		{
			nx = eigVecs.ptr<double>(0)[2];
			ny = eigVecs.ptr<double>(0)[3];
			eigValue = eigs.ptr<double>(0)[1];
		}

		double temp = sqrt(nx * nx + ny * ny);
		if (temp != 0.0)
		{
			// 法向量先进行单位化
			nx = nx / temp;
			ny = ny / temp;
			double t = -(rx * nx + ry * ny) / (rxx * nx * nx + 2 * rxy * nx * ny + ryy * ny * ny);
			double px = t * nx;
			double py = t * ny;
			if (fabs(px) <= 0.5 && fabs(py) <= 0.5)
			{
				idTmp[m] = m;
				dTmp[m] = cv::Point2f(j - px, i - py);
				//d2.push_back(cv::Point2f(j - px, i - py));
			}
		}
	}

	for (int i = 0; i < count; i++)
	{
		if (idTmp[i] != -1)
		{
			d2.push_back(dTmp[i]);
		}
	}


	cv::RotatedRect minEllipse = fitEllipse(cv::Mat(pts));
	centerPt.x = minEllipse.center.x;
	centerPt.y = minEllipse.center.y;

	// 要是点太少 就木有亚像素的必要了，直接返回
	if (d2.size() < m_minNumpt || d2.size() < 6) {
		return false;
	}

	// 去除椭圆边沿上误差较大的点
	minEllipse = fitEllipse(cv::Mat(d2));
	double centerX = minEllipse.center.x;
	double centerY = minEllipse.center.y;
	double c = sqrt(fabs(minEllipse.size.width * minEllipse.size.width - minEllipse.size.height * minEllipse.size.height)) / 2.0;
	//printf("center before = %.7f %.7f\n", centerX, centerY);

	//double angle = AngleAxes(minEllipse);

	double angle = minEllipse.angle;
	if (angle > 180) {
		angle = angle - 180;
	}
	if (minEllipse.size.width < minEllipse.size.height)
		angle = angle - 90;

	angle = angle * CV_PI / 180.0;

	double leftCenter_x2 = centerX - c * cos(angle);
	double leftCenter_y2 = centerY - c * sin(angle);
	double rightCenter_x2 = centerX + c * cos(angle);
	double rightCenter_y2 = centerY + c * sin(angle);

	vector<cv::Point2f> d3;
	for (int i = 0; i < (int)d2.size();)
	{
		double ptError = fabs((distance(leftCenter_x2, leftCenter_y2, d2[i].x, d2[i].y)
			+ distance(rightCenter_x2, rightCenter_y2, d2[i].x, d2[i].y) - minEllipse.size.height));

		if (ptError > FitErrorThreshold) {
			i++;
			continue;
		}
		d3.push_back(d2[i]);
		i++;
	}

	// 筛完之后统计点的个数 个数少就不亚像素啦！
	if (d3.size() < m_minNumpt || d3.size() < 6) {
		return false;
	}

	minEllipse = fitEllipse(cv::Mat(d3));
	centerPt.x = minEllipse.center.x;
	centerPt.y = minEllipse.center.y;
	//printf("center after = %.7f %.7f\n", centerPt.x, centerPt.y);
	//cout << "centerPt after = " << centerPt << endl;
	return true;
}


// 两点间【距离】
double distance(double x1, double y1, double x2, double y2) {
	double dx = x2 - x1;
	double dy = y2 - y1;
	return (sqrt(dx * dx + dy * dy));
}