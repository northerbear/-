if __name__ == '__main__':
img_path_3 = 'yuelian.png'
// ��ȡͼ��
image = cv2.imread(img_path_3, cv2.IMREAD_GRAYSCALE)
image = cv2.resize(image, (512, 512))

//����x��y�����ϵ��ݶ�
sobel_x = cv2.Sobel(image, cv2.CV_64F, 1, 0, ksize = 3)
sobel_x = cv2.convertScaleAbs(sobel_x)
sobel_y = cv2.Sobel(image, cv2.CV_64F, 0, 1, ksize = 3)
sobel_x = cv2.convertScaleAbs(sobel_y)
img_show('==', image)
img_show('==', sobel_x)
img_show('==', sobel_y)
// ���ݶ�ͼ��ϲ�Ϊһ���ۺϵı�Եͼ��
sobel_combined = cv2.addWeighted(sobel_x.astype(np.uint8), 0.5, sobel_y.astype(np.uint8), 0.5, 0)
#//��ʾԭʼͼ��ͱ�Ե�����
img_show('--', sobel_combined)

