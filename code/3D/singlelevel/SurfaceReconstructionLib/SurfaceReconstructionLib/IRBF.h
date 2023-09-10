#pragma once
#include "Vec3d.h"
#include "Mat3.h"
#include "Octree.h"
#include <vector>
#include <cmath>
#include <time.h>
#include <float.h>
#include <queue>
#include"polygonizer.h"
#include "nearest.h"
#include <iostream>

#define EPS 1e-7


class IRBF : public ImplicitFunction {
public:
	IRBF(double r, double s, int n, double mu,Vec3d* points, Vec3d* normals, Octree* octree,double c);
	~IRBF();
	inline virtual double eval(double x, double y, double z);
	void SurfaceReconstruction(int res, Vec3d minPoint, Vec3d maxPoint);  // �����ؽ�
	double* resultValues;  // �����ؽ����((res+1)^3)
private:
	int maxk;  // �����ڸ���
	double r;//ƽ��֧�Ű뾶
	double s;  // CSRBF����
	double mu;
	int res;  // �����ؽ�����
	Vec3d minPoint, maxPoint;  // �ؽ���Χ
	double c;

	int n;  // �����ڵ����
	Vec3d* points;  // ������������(�����������ⲿ����)
	Vec3d* normals;  // ���Ʒ���������(�����������ⲿ����)
	Octree* octree;  // ����(�����������ⲿ����)

	Mat3* orthbases;  // ���Ƹ��㴦�ֲ�����ϵ�Ļ�
	Vec3d* hparas;  // ���Ƹ����Ӧ�����������
	double* lserror;  // ���Ƹ������������ϵ���С�������
	double* approxArgs;  // ���Ƹ����Ӧ�����ϵ��
	double* cis;  // ���Ƹ����Ӧ��Ci
	double boxlength;


	Nearest* nearest; //������Ϣ����СΪn
	double* lambda;
	double* nearestdistance;
	double* mnearestdistance;
	int averagecount;

	inline double RBF(const double& squaredDistanceRatio, const double& approxArg) const;  // ���ϵ�����Ƶĵ�RBF����Ȼ�����
	inline void CSRBF(const std::vector<int>* Indices, const std::vector<double>* squaredDists, double* values,int k) const;  // ����������볡(��С����)��CSRBF����ֵ��
	inline double CalHvalue(const Vec3d& point, const int& i) const;  // ����������ڵ��Ƶ�i�����Ӧ�ľֲ����������µĶ�Ӧhֵ

	void FindNearestNeighborsAndCalculateHparas(); // ������ƵĽ�����Ϣ��ÿ����Ķ����������
	void CalculateApproxArgs();  // �������������ϵ��
	void CalculateCis();  // Ϊ���Ƽ���ÿ�����Ciֵ
	void CalculateReconstructionValues();  // ���������ؽ����
	void DeleteNearest();
	void CaculateInditor();

	
};



void IRBF::DeleteNearest()
{
	for(int i=0;i<n;i++)
		nearest[i].Deletenearest();
}

IRBF::IRBF(double r, double s, int n, double mu,Vec3d* points, Vec3d* normals, Octree* octree,double c) {
	// �������
	this->r = r;
	this->s = s;
	this->n = n;
	this->mu = mu;
	this->c = c;
	this->points = points;
	this->normals = normals;
	this->octree = octree;
	orthbases = new Mat3[n];
	hparas = new Vec3d[n];
	lserror = new double[n];
	approxArgs = new double[n];
	cis = new double[n];
	nearestdistance = new double[n];
	mnearestdistance = new double[n];
	lambda = new double[n];
	maxk = 0;
	averagecount = 0;
	boxlength = (maxPoint - minPoint).SquaredLength();
	// ��ʼ�������ؽ��������Ϣ
	printf("Precalculating info for surface reconstruction...\n");
	clock_t start = clock();

	nearest = new Nearest[n];
	FindNearestNeighborsAndCalculateHparas();
	CaculateInditor();
	CalculateApproxArgs();
	CalculateCis();



	double cost = (double)(clock() - start) / 1000;
	printf("Finished in %.3lf seconds.\n", cost);

	DeleteNearest();


}

IRBF::~IRBF() {
	delete[] orthbases;
	delete[] hparas;
	delete[] lserror;
	delete[] approxArgs;
	delete[] cis;
	delete[] nearestdistance;
	delete[] mnearestdistance;
	delete[] lambda;
}


void IRBF::SurfaceReconstruction(int res, Vec3d minPoint, Vec3d maxPoint) {
	this->res = res;
	this->minPoint = minPoint;
	this->maxPoint = maxPoint;
	if (resultValues != NULL) { delete[] resultValues; }
	resultValues = new double[(res + 1) * (res + 1) * (res + 1)];
	CalculateReconstructionValues();
}


inline double IRBF::RBF(const double& rho, const double& approxArg) const {
	double sqrtr = sqrt(rho); double temp = (1 - sqrtr) * (1 - sqrtr);
	return temp * temp * (4 * sqrtr + 1) * (1 / sqrt(approxArg + rho) );
}

inline void IRBF::CSRBF(const std::vector<int>* Indices, const std::vector<double>* squaredDists, double* values,int k) const {
	if ((*squaredDists)[0] < EPS && approxArgs[(*Indices)[0]] < EPS) {
		values[0] = 1;
		for (int i = 1; i < k; i++) {
			values[i] = 0;
		}
		return;
	}
	else
	{
		double sumVal = 0;
		for (int i = 0; i < k; i++) {
			values[i] = RBF(lambda[(*Indices)[i]]*(*squaredDists)[i] / r, approxArgs[(*Indices)[i]]);
			sumVal += values[i];
		}
		for (int i = 0; i < k; i++) {
			values[i] /= sumVal;
		}
	}

}

inline double IRBF::CalHvalue(const Vec3d& point, const int& i) const {
	Vec3d transformedLoc = orthbases[i].BaseTransform(point - points[i]);
	double surfaceW =
		hparas[i].x * transformedLoc.x * transformedLoc.x +
		hparas[i].y * transformedLoc.x * transformedLoc.y +
		hparas[i].z * transformedLoc.y * transformedLoc.y;
	return (transformedLoc.z - surfaceW);
}


inline void IRBF::FindNearestNeighborsAndCalculateHparas() {
 #pragma omp parallel for
	for (int i = 0; i < n; i++) {
		Nearest t_nearest = nearest[i];
		t_nearest.count = octree->FindNearestNeighborsBySquaredRadius(points[i], r, t_nearest.indices, t_nearest.squaredDists, nearestdistance[i]);
		// ����õ㴦�ֲ�����ϵ�Ļ�
		mnearestdistance[i] = ((*t_nearest.squaredDists)[1]);
		nearestdistance[i] /= r;
		//lambda[i] = sqrt(mu / mnearestdistance[i]);
		//std::cout << lambda[i] << std::endl;
		lambda[i] = 1;
		nearest[i].count = t_nearest.count;
		Vec3d tempvec = normals[i];
		tempvec.ApplyAbs();
		if ((tempvec.x <= tempvec.y) && (tempvec.x <= tempvec.z)) {
			tempvec = Vec3d(1, 0, 0);
		}
		else if ((tempvec.y <= tempvec.x) && (tempvec.y <= tempvec.z)) {
			tempvec = Vec3d(0, 1, 0);
		}
		else {
			tempvec = Vec3d(0, 0, 1);
		}
		Vec3d base1 = normals[i].Cross(tempvec);
		Vec3d base2 = normals[i].Cross(base1);
		orthbases[i] = Mat3(base1, base2, normals[i]);
		// ʹ����С���˷�����õ㴦��������ϵ��
		Mat3 tempmat = Mat3();
		tempvec = Vec3d();
		for (int j = 0; j < t_nearest.count; j++) {
			Vec3d t;
			t = points[(*t_nearest.indices)[j]] - points[i];
			Vec3d transformedLoc = orthbases[i].BaseTransform(t);
			tempmat.u.x += transformedLoc.x*transformedLoc.x*transformedLoc.x*transformedLoc.x;
			tempmat.v.x += transformedLoc.x*transformedLoc.x*transformedLoc.x*transformedLoc.y;
			tempmat.w.x += transformedLoc.x*transformedLoc.x*transformedLoc.y*transformedLoc.y;
			tempmat.w.y += transformedLoc.x*transformedLoc.y*transformedLoc.y*transformedLoc.y;
			tempmat.w.z += transformedLoc.y*transformedLoc.y*transformedLoc.y*transformedLoc.y;
			tempvec.x += transformedLoc.x*transformedLoc.x*transformedLoc.z;
			tempvec.y += transformedLoc.x*transformedLoc.y*transformedLoc.z;
			tempvec.z += transformedLoc.y*transformedLoc.y*transformedLoc.z;
		}
		tempmat.u.y = tempmat.v.x;
		tempmat.u.z = tempmat.v.y = tempmat.w.x;
		tempmat.v.z = tempmat.w.y;
		hparas[i] = tempmat.SolveLinerEquation(tempvec);  // ����������ϵ��

	}
	//std::cout << minlserror << '\t' << maxlserror << std::endl;
}


inline void IRBF::CaculateInditor() 
{
#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		Nearest t_nearest = nearest[i];
		t_nearest.count = octree->FindNearestNeighborsBySquaredRadius(points[i], r, t_nearest.indices, t_nearest.squaredDists, nearestdistance[i]);
		mnearestdistance[i] = sqrt((*t_nearest.squaredDists)[1]);
		nearestdistance[i] /= r;
		nearest[i].count = t_nearest.count;
		//lambda[i] = 1;
		maxk = t_nearest.count > maxk ? t_nearest.count : maxk;


		// ������С���˷�������
		double err = 0;
		for (int j = 0; j < t_nearest.count; j++) {
			int j1 = (*t_nearest.indices)[j];
			double temp = CalHvalue(points[i], j1);
			double r1 = sqrt((*t_nearest.squaredDists)[j]/nearest[j1].r); double r2 = (1 - r1) * (1 - r1);
			err += temp*temp*r2*r2;
		}
		lserror[i] = sqrt(err);
	}
}


void IRBF::CalculateApproxArgs() {
	double meanerr = 0;double maxerr = 0, minerr = 0;
	int i, j;
	//std::cout << maxlserror << '\t' << minlserror << std::endl;
	// ��mink��max
	for (i = 0; i < n; i++) {
		meanerr += lserror[i];
		maxerr = maxerr > lserror[i] ? maxerr : lserror[i];
		minerr = minerr < lserror[i] ? minerr : lserror[i];
	}
	meanerr = meanerr / n; 
	//std::cout << meanerr<<varerr << std::endl;

	std::cout << "center value:"<<c << std::endl;
#pragma omp parallel for
	for (j = 0; j < n; j++) {
		approxArgs[j] = c;
		if (lserror[j] - meanerr > EPS)
		{
			double temp = (lserror[j] - meanerr) / (maxerr - minerr);
			approxArgs[j] = c/(1+temp);
		}
		else if (lserror[j] - meanerr < -EPS)
		{
			double temp = (lserror[j] - meanerr) / (maxerr - minerr);
			approxArgs[j] = c * (1-temp);
		}
		//approxArgs[j] = c;
	}
}

void IRBF::CalculateCis() {
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		double* CSRBFValues = new double[maxk];  // �ݴ�CSRBF����ֵ��
		Nearest t_nearest = nearest[i];
		CSRBF(t_nearest.indices, t_nearest.squaredDists, CSRBFValues, t_nearest.count);
		double c = 0;
		for (int j = 0; j < t_nearest.count; j++) {
			c -= CalHvalue(points[i], (*t_nearest.indices)[j]) * CSRBFValues[j];
		}
		cis[i] = c;
		delete[] CSRBFValues;
	}

}

inline double IRBF::eval(double x, double y, double z) {
	Vec3d point(x, y, z);
	std::vector<int>* indices = new std::vector<int>;
	std::vector<double>* squaredDists = new std::vector<double>;
	int findCount = octree->FindNearestNeighborsBySquaredRadius(point, r, indices, squaredDists);
	double chvsum = 0, vsum = 0, rbfValue; double t_approx;
	int legalCount = 0;  // ʶ��õ����������ݵ��������
	for (int i = 0; i < findCount; i++) {
		if ((*squaredDists)[i] >= r) {
			continue;  // �õ㲻�����ݵ�֧�Ű뾶��������
		}
		legalCount++;
		t_approx = approxArgs[(*indices)[i]];
		if (t_approx < EPS && (*squaredDists)[i] < EPS) {
			//approxArgs[(*indices)[i]<EPS��ʾ�õ�Ϊ���Ƶ�
			chvsum = 0;
			vsum = 1;
			break;
		}
		rbfValue = RBF( lambda[(*indices)[i]]*(*squaredDists)[i] / r, t_approx);
		vsum += rbfValue;
		chvsum += (cis[(*indices)[i]] + CalHvalue(point, (*indices)[i])) * rbfValue;
	}
	if (legalCount < 1) {  // �õ�õ����������ݵ������������
		return(32768);
	}
	else {
		return(chvsum / vsum);
	}
	//return x * x + y * y + z * z - 1;
	delete indices;
	delete squaredDists;

}

void IRBF::CalculateReconstructionValues() {
	printf("|Calculating--------------------------------------|\n");
	// ����������
	int progressbar_step = res / 50;
	// ��ʱ
	clock_t start = clock();
	Vec3d stepVec = (maxPoint - minPoint) / res;
	int count = res + 1;
	int sqcount = count * count;
#pragma omp parallel for
	for (int iz = 0; iz <= res; iz++) {
		std::vector<int>* indices = new std::vector<int>();
		std::vector<double>* squaredDists = new std::vector<double>();
		for (int ix = 0; ix <= res; ix++) {
			for (int iy = 0; iy <= res; iy++) {
				Vec3d point(ix, iy, iz);
				point.ApplyMultiplication(stepVec);
				point += minPoint;
				int findCount = octree->FindNearestNeighborsBySquaredRadius(point, r, indices, squaredDists);
				int index = iz * sqcount + ix * count + iy;
				// ����õ����Ϻ���ֵ
				double chvsum = 0, vsum = 0, rbfValue;
				int legalCount = 0;  // ʶ��õ����������ݵ��������
				for (int i = 0; i < findCount; i++) {
					if ((*squaredDists)[i] >= r) {
						continue;  // �õ㲻�����ݵ�֧�Ű뾶��������
					}
					legalCount++;
					if (approxArgs[(*indices)[i]] < EPS && (*squaredDists)[i] < EPS) {
						chvsum = 0;
						vsum = 1;
						break;
					}
					rbfValue = RBF(lambda[(*indices)[i]]*(*squaredDists)[i] / r, approxArgs[(*indices)[i]]);
					vsum += rbfValue;
					chvsum += (cis[(*indices)[i]] + CalHvalue(point, (*indices)[i])) * rbfValue;
				}
				if (legalCount < 1) {  // �õ�õ����������ݵ������������
					resultValues[index] = 32768;
				}
				else {
					resultValues[index] = chvsum / vsum;
				}

			}
		}
		delete indices;
		delete squaredDists;
		if (iz % progressbar_step == 0) { printf(">"); }
	}

	double cost = (double)(clock() - start) / 1000;
	printf("\nFinished in %.3lf seconds.\n", cost);
}

