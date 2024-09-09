if __name__ == '__main__':
img_path_3 = 'yuelian.png'
// 读取图像
image = cv2.imread(img_path_3, cv2.IMREAD_GRAYSCALE)
image = cv2.resize(image, (512, 512))

//计算x和y方向上的梯度
sobel_x = cv2.Sobel(image, cv2.CV_64F, 1, 0, ksize = 3)
sobel_x = cv2.convertScaleAbs(sobel_x)
sobel_y = cv2.Sobel(image, cv2.CV_64F, 0, 1, ksize = 3)
sobel_x = cv2.convertScaleAbs(sobel_y)
img_show('==', image)
img_show('==', sobel_x)
img_show('==', sobel_y)
// 将梯度图像合并为一个综合的边缘图像
sobel_combined = cv2.addWeighted(sobel_x.astype(np.uint8), 0.5, sobel_y.astype(np.uint8), 0.5, 0)
#//显示原始图像和边缘检测结果
img_show('--', sobel_combined)

