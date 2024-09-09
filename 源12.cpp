




std::cout << "��֪��Դ��\n" << emarkerPoints << std::endl;

Eigen::MatrixXd x(3, 3);
x << 1, 0, 0,
0, 0.9396926, -0.3420201,
0, 0.3420201, 0.9396926;
Eigen::MatrixXd y(3, 3);
y << 0.8660254, 0, 0.5,
0, 1, 0,
-0.5, 0, 0.8660254;
Eigen::MatrixXd z(3, 3);
z << 0.7660444, -0.6427876, 0,
0.6427876, 0.7660444, 0,
0, 0, 1;
Eigen::MatrixXd Rr(3, 3);
Rr = x * y * z;
std::cout << "������ת����\n" << Rr << std::endl;


Eigen::MatrixXd etransPoints(lPoints.size(), 3);  //����ת����ĵ��������
Eigen::MatrixXd transform(3, 4);  //rt�ϲ���ת������
transform << 0.66341391, -0.55667039, 0.5, 9,
0.73502404, 0.60992311, -0.29619809, 5,
-0.14007685, 0.56401396, 0.81379766, 6;
Eigen::MatrixXd eesourcePoints(lPoints.size(), 4); // Դ������Ϊ���������ʽ

int i, j;
for (i = 0; i < lPoints.size(); i++)
{
	eesourcePoints.row(i) << emarkerPoints.row(i), 1;
}
etransPoints.transpose() = transform * eesourcePoints.transpose();//�õ��ڶ����ʼ����

//double sigma = 0.05;
//addGaussianNoise(etransPoints, sigma);
//���ϸ�˹����


srand(static_cast<unsigned int>(time(nullptr))); //c��ʼ��
for (i = 0; i < lPoints.size(); i++) {
	for (j = 0; j < 3; j++) {
		int r = rand();
		double gaosi = static_cast<double>(r) / RAND_MAX * 0.5;
		etransPoints(i, j) += gaosi;
	}
}

//ת����ĵ����
std::cout << "����ڶ�����ά��Ϊ��\n" << etransPoints << std::endl;


Eigen::Matrix3d R;
Eigen::Vector3d t;
svdRegistration(emarkerPoints, etransPoints, R, t);

std::cout << "Rotation Matrix: \n" << R << std::endl;
std::cout << "Translation Vector: \n" << t << std::endl;
// �õ�rt
Eigen::MatrixXd T(3, 4);
for (i = 0; i < 3; i++)
{
	T.row(i) << R.row(i), t.row(i);
}

Eigen::MatrixXd erorrPoints(3, 4);
erorrPoints = transform - T;
std:cout << "����С��\n" << erorrPoints << std::endl;




return 0;









std::cout << "����ת����Դ�����CΪ��\n" << emarkerPoints << std::endl;
std::cout << "���ڻ�׼�ľ���BΪ��\n" << esourcePoints << std::endl;



ifstream fileDst, filedst2;
vector<glm::dvec3> dst, dstr;
Eigen::MatrixXd ePoints(64, 3);
Eigen::MatrixXd mPoints(64, 3);

int i, j;
filedst2.open("c.real.txt");
assert(filedst2.is_open());   //��ʧ��,�����������Ϣ,����ֹ�������� 
while (!filedst2.eof())
{
	for (i = 0; i < 64; i++)
	{
		for (j = 0; j < 3; j++) {
			filedst2 >> mPoints(i, j);
		}
	}
}

std::cout << "ʵ��ת���㣺\n" << mPoints << std::endl;

fileDst.open("B.real.txt");   //���ļ����������ļ��������� 
assert(fileDst.is_open());   //��ʧ��,�����������Ϣ,����ֹ�������� 
while (!fileDst.eof())
{
	for (i = 0; i < 64; i++)
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

Eigen::MatrixXd transform(3, 4);  //rt�ϲ���ת������

for (i = 0; i < 3; i++)
{
	transform.row(i) << R.row(i), t.row(i);
}

Eigen::MatrixXd etransPoints(64, 3);  //����ת����ĵ��������
Eigen::MatrixXd eesourcePoints(64, 4); // Դ������Ϊ���������ʽ
//Eigen::MatrixXd ematchPoints(70, 4); // Դ������Ϊ���������ʽ
for (i = 0; i < 64; i++)
{
	eesourcePoints.row(i) << mPoints.row(i), 1;
	//ematchPoints.row(i) << mPoints.row(i), 1;
}

etransPoints.transpose() = transform * eesourcePoints.transpose();   //ת����ĵ����
//etransPoints.transpose() = transform * ematchPoints.transpose();

std::cout << "ͨ��ת������õ�����ά��Ϊ��\n" << etransPoints << std::endl;


Eigen::MatrixXd eorrorPoints(64, 3);  //����������
//eorrorPoints = mPoints - etransPoints;
eorrorPoints = ePoints - etransPoints;

std::cout << "The error is :\n" << eorrorPoints << std::endl;



return 0;

















Eigen::MatrixXd transform(3, 4);  //rt�ϲ���ת������

for (i = 0; i < 3; i++)
{
	transform.row(i) << R.row(i), t.row(i);
}

Eigen::MatrixXd etransPoints(64, 3);  //����ת����ĵ��������
Eigen::MatrixXd eesourcePoints(64, 4); // Դ������Ϊ���������ʽ
//Eigen::MatrixXd ematchPoints(70, 4); // Դ������Ϊ���������ʽ
for (i = 0; i < 64; i++)
{
	eesourcePoints.row(i) << mPoints.row(i), 1;
	//ematchPoints.row(i) << mPoints.row(i), 1;
}

etransPoints.transpose() = transform * eesourcePoints.transpose();   //ת����ĵ����
//etransPoints.transpose() = transform * ematchPoints.transpose();

std::cout << "ͨ��ת������õ�����ά��Ϊ��\n" << etransPoints << std::endl;


Eigen::MatrixXd eorrorPoints(64, 3);  //����������
//eorrorPoints = mPoints - etransPoints;
eorrorPoints = ePoints - etransPoints;

std::cout << "The error is :\n" << eorrorPoints << std::endl;



return 0;

}
