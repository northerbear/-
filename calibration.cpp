// 双目相机标定及校正.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include<opencv2/core/core.hpp>
#include<opencv2/calib3d/calib3d.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<iostream>
#include<fstream>
#define CV_RGB2GRAY cv::COLOR_RGB2GRAY

using namespace std;
using namespace cv;



void realPointsCreate(vector<String> imgPathList, vector<vector<Point3f>>&Points_object);//获取标定板特征点在世界坐标系下的三维坐标
void draw(Mat img1, Mat img2, Size imgsize, Rect validroi[2]);//画出等距极限检验图像校正结果
bool singleCameraCalibrate(vector<String> imgPathList, const char* singleCalibrate_result, vector<vector<Point3f>> Points_object,
    vector<vector<Point2f>>& Points, Mat& cameraK, Mat& distCoeffs, vector<Mat>& rotationMat, vector<Mat>& translationMat, 
    Size& imagesize, Size& board_size, Size& block_size);//单目相机标定并保存标定结果到指定文件



vector<String>imgPathList_L;//用于标定的图像列表
vector<String>imgPathList_R;
const char* image_temp_L = "E://image//深度测量//image7-25//calibration//left//10000_0.bmp";//用于检测深度的图像
const char* image_temp_R = "E://image//深度测量//image7-25//calibration//right//10000_1.bmp";
const char* singleCalibrate_result_L = "calibrationresult_L.txt";//保存单目标定的结果
const char* singleCalibrate_result_R = "calibrationresult_R.txt";
const char* stereoCalibrate_result = "stereoCalibrateResult.txt";//保存立体标定的结果
const char* stereoRectifyParams = "stereoRectifyParams.txt";//保存立体校正参数
vector<vector<Point2f>>points_L;//保存标定板特征点在相机像素坐标系下的二维坐标
vector<vector<Point2f>>points_R;
vector<vector<Point3f>>points_object;//标定板特征点在世界坐标系的三维坐标
Mat cameraK_L = Mat(3, 3, CV_32FC1, Scalar::all(0)); //相机内参矩阵3*3
Mat cameraK_R = Mat(3, 3, CV_32FC1, Scalar::all(0));
Mat distCoeffs_L = Mat(1, 5, CV_32FC1, Scalar::all(0)); //相机畸变矩阵1*5
Mat distCoeffs_R = Mat(1, 5, CV_32FC1, Scalar::all(0));
vector<Mat> rotation_L; //旋转矩阵 
vector<Mat> rotation_R; 
vector<Mat> translation_L; //平移矩阵
vector<Mat> translation_R; 
Mat R, T, E, F;//立体标定参数
Mat R1, R2, P1, P2, Q;//立体校正参数
Mat maplx, maply, maprx, mapry;//图像重投影映射表
Mat imgL_rectified, imgR_rectified, disparity, img_result3D;//校正图像 视差图 深度图
Size board_size = Size(9, 6);//标定板每行每列角点个数
Size block_size = Size(25, 25);  //标定板上每个方格的实际大小，只影响最后求解的平移向量
Size imagesize;//图像尺寸
Rect validroi[2];//校正图像裁剪后的区域



int main()
{
    //左右相机分别进行标定
    glob("E://image//深度测量//image7-25//calibration//left", imgPathList_L);
    glob("E://image//深度测量//image7-25//calibration//right", imgPathList_R);
    realPointsCreate(imgPathList_L, points_object);
    cout << "开始对左相机进行标定" << endl;
    singleCameraCalibrate(imgPathList_L, singleCalibrate_result_L, points_object, points_L, cameraK_L, distCoeffs_L, rotation_L, translation_L,
        imagesize, board_size, block_size);
    cout << "开始对右相机进行标定" << endl;
    singleCameraCalibrate(imgPathList_R, singleCalibrate_result_R, points_object, points_R, cameraK_R, distCoeffs_R, rotation_R, translation_R,
        imagesize, board_size, block_size);


    //立体标定
    cout << "开始立体标定并保存标定结果---------->";
    ofstream stereo;
    stereo.open(stereoCalibrate_result);
    TermCriteria criteria = TermCriteria(TermCriteria::COUNT + TermCriteria::EPS, 30, 1e-6);//迭代终止条件
    stereoCalibrate(points_object, points_L, points_R, cameraK_L, distCoeffs_L, cameraK_R, distCoeffs_R, imagesize, 
        R, T, E, F, CALIB_FIX_INTRINSIC, criteria);
    stereo << "旋转矩阵" << endl;
    stereo << R << endl;
    stereo << "平移矩阵" << endl;
    stereo << T << endl;
    stereo << "本质矩阵" << endl;
    stereo << E << endl;
    stereo << "基础矩阵" << endl;
    stereo << F << endl;
    stereo.close();
    cout << "success" << endl;


    //立体校正
    cout << "创建图像重投影映射表---------->";
    stereoRectify(cameraK_L, distCoeffs_L, cameraK_R, distCoeffs_R, imagesize, R, T, R1,
        R2, P1, P2, Q, 0, -1, imagesize, &validroi[0], &validroi[1]);
    stereo.open(stereoRectifyParams);
    stereo << "R1" << endl;
    stereo << R1 << endl;
    stereo << "R2" << endl;
    stereo << R2 << endl;
    stereo << "P1" << endl;
    stereo << P1 << endl;
    stereo << "P2" << endl;
    stereo << P2 << endl;
    stereo << "Q" << endl;
    stereo << Q << endl;
    stereo.close();
    
    Mat img_L = imread(image_temp_L);
    Mat img_R = imread(image_temp_R);
    if (img_L.empty()||img_R.empty())
    {
        cout << "错误；请打开正确的图片";
        return -1;
    }
   
    initUndistortRectifyMap(cameraK_L, distCoeffs_L, R1, P1, imagesize, CV_32FC1, maplx, maply);
    initUndistortRectifyMap(cameraK_R, distCoeffs_R, R2, P2, imagesize, CV_32FC1, maprx, mapry);
    cout << "success" << endl;
    cout << "图像校正并画出等距极线进行检验" << endl;
    remap(img_L, imgL_rectified, maplx, maply, INTER_LINEAR);
    remap(img_R, imgR_rectified, maprx, mapry, INTER_LINEAR);
    draw(imgL_rectified, imgR_rectified, imagesize, &validroi[2]);

    //立体匹配
    cv::cvtColor(imgL_rectified, imgL_rectified, CV_RGB2GRAY);
    cv::cvtColor(imgR_rectified, imgR_rectified, CV_RGB2GRAY);
    cout << "创建视差图" << endl;
    Ptr<StereoBM> bm = StereoBM::create(16, 9);//智能指针
    bm->compute(imgL_rectified, imgR_rectified, disparity);//计算视差图
    disparity.convertTo(disparity, CV_32F, 1.0 / 16);
    normalize(disparity, disparity, 0, 256, NORM_MINMAX, -1);//归一化视差映射
    reprojectImageTo3D(disparity, img_result3D, Q); 
    cout << "完成" << endl;
    imshow("视差图", disparity);
    imshow("深度图", img_result3D);
    //得到的深度图为三通道，每个像素具有三个通道，分别每个像素的在相机坐标系下的三维坐标。
  
    waitKey(0);
    return 0;
}


void draw(Mat imgL_rectified, Mat imgR_rectified, Size imagesize, Rect validroi[2])
{
    cv::Mat canvas;
    double sf;
    int w, h;
    sf = 600. / MAX(imagesize.width, imagesize.height);
    w = cvRound(imagesize.width * sf);
    h = cvRound(imagesize.height * sf);
    canvas.create(h, w * 2, CV_8UC3);

    /*左图像画到画布上*/
    Mat canvasPart = canvas(Rect(w * 0, 0, w, h));
    resize(imgL_rectified, canvasPart, canvasPart.size(), 0, 0, INTER_AREA);//把图像缩放到跟canvasPart一样大小  
    Rect vroiL(cvRound(validroi[0].x * sf), cvRound(validroi[0].y * sf),//获得被截取的区域    
        cvRound(validroi[0].width * sf), cvRound(validroi[0].height * sf));
    rectangle(canvasPart, vroiL, Scalar(0, 0, 255), 3, 8); //画上一个矩形  

    /*右图像画到画布上*/
    canvasPart = canvas(Rect(w, 0, w, h)); //获得画布的另一部分  
    resize(imgR_rectified, canvasPart, canvasPart.size(), 0, 0, INTER_LINEAR);
    Rect vroiR(cvRound(validroi[1].x * sf), cvRound(validroi[1].y * sf),
        cvRound(validroi[1].width * sf), cvRound(validroi[1].height * sf));
    rectangle(canvasPart, vroiR, Scalar(0, 255, 0), 3, 8);

    /*画上对应的线条*/
    for (int i = 0; i < canvas.rows; i += 16)
        line(canvas, Point(0, i), Point(canvas.cols, i), Scalar(0, 255, 0), 1, 8);
    imwrite("result.bmp", canvas);
    imshow("rectified", canvas);
}


void realPointsCreate(vector<String> imgPathList, vector<vector<Point3f>>&Points_object)
{
    for (int t = 0; t < imgPathList.size(); t++)
    {
        vector<Point3f> points3D_per_img;
        for (int i = 0; i < board_size.height; i++)
        {
            for (int j = 0; j < board_size.width; j++)
            {
                //假定标定板放在世界坐标系中z=0的平面上
                Point3f realPoints;
                realPoints.x = i * block_size.width;
                realPoints.y = j * block_size.height;
                realPoints.z = 0;
                points3D_per_img.push_back(realPoints);
            }
        }
        Points_object.push_back(points3D_per_img);
    }//初始化标定板上角点的三维坐标
}



//单目标定
bool singleCameraCalibrate(vector<String> imgPathList, const char* singleCalibrate_result, vector<vector<Point3f>> Points_object,
    vector<vector<Point2f>>& Points, Mat& cameraK, Mat& distCoeffs, vector<Mat> &rotationMat, vector<Mat> &translationMat,
    Size& imagesize, Size& board_size, Size& block_size)
{
    int img_nums = 0;//用于标定的图片的数量
    vector<Point2f>points_per_img;//存放每幅图像检测到的角点
    for (img_nums = 0; img_nums < imgPathList.size(); img_nums++)
    {
        Mat img = imread(imgPathList[img_nums]);
        Mat dst;
        cvtColor(img, dst, CV_RGB2GRAY);

        if (img_nums == 0)   //读入第一张图片时 获取图片的宽高信息
        {
            imagesize.height = img.rows;
            imagesize.width = img.cols;
            cout << "图片的宽度为：" << imagesize.width << endl;
            cout << "图片的高度为：" << imagesize.height << endl;
            cout << endl;
            cout << "开始提取角点" << endl;
        }


        //step1：提取角点
        cout << "正在提取第" << img_nums + 1 << "幅图的角点---------->";
        bool success = findChessboardCorners(dst, board_size, points_per_img);
        if (!success)
        {
            cout << "can not find chessboard corners" << endl;
            return false;
        }
        else
        {
            //step2亚像素精确化（两种方法）
            find4QuadCornerSubpix(dst, points_per_img, Size(5, 5));
            //cornerSubPix(img, points_per_img, Size(5, 5));
            Points.push_back(points_per_img);//保存亚像素角点
            //step3在图像上显示角点位置
            drawChessboardCorners(img, board_size, points_per_img, success);
            cout << "success" << endl;
        }

    }cout << "角点提取完成" << endl << endl;


    //step4 标定
    cout << "开始标定---------->";
    int point_counts = board_size.area();//初始化每个图像中角点的数量
    calibrateCamera(Points_object, Points, imagesize, cameraK, distCoeffs, rotationMat, translationMat, 0);
    cout << "标定完成" << endl;

    //step5评定标定结果
    cout << "开始评价标定结果" << endl;
    double total_err = 0.0;//所有图像的平均误差
    double err = 0.0;//每个图像的平均误差
    vector<Point2f>repoints_per_img;//重新保存计算得到的新的投影点
    for (int i = 0; i < img_nums; i++)
    {
        //通过得到的摄像机内外参数，对空间的三维点进行重新投影计算，得到新的投影点
        projectPoints(Points_object[i], rotationMat[i], translationMat[i], cameraK, distCoeffs, repoints_per_img);
        //计算新的投影点和旧的投影点之间的误差
        vector<Point2f> tempImagePoint = Points[i];
        Mat tempImagePointMat = Mat(1, tempImagePoint.size(), CV_32FC2);
        Mat image_points2Mat = Mat(1, repoints_per_img.size(), CV_32FC2);
        for (int j = 0; j < tempImagePoint.size(); j++)
        {
            image_points2Mat.at<Vec2f>(0, j) = Vec2f(repoints_per_img[j].x, repoints_per_img[j].y);
            tempImagePointMat.at<Vec2f>(0, j) = Vec2f(tempImagePoint[j].x, tempImagePoint[j].y);
        }
        err = norm(image_points2Mat, tempImagePointMat, NORM_L2);
        total_err += err /= point_counts;//54
        cout << "第" << i + 1 << "幅图像的平均误差：" << err << "像素" << endl;
    }cout << "平均误差：" << total_err / img_nums << "像素" << endl << "评价完成" << endl;//17

    //开始保存相机的参数
    ofstream result(singleCalibrate_result);
    cout << "开始保存相机的标定结果" << endl;
    result << "相机的内参矩阵" << endl;
    result << cameraK << endl << endl;
    result << "相机的畸变参数" << endl;
    result << distCoeffs << endl;
    cout << "保存成功" << endl;
    cout << "标定完成" << endl << endl << endl;
    result.close();
    
    return true;
}