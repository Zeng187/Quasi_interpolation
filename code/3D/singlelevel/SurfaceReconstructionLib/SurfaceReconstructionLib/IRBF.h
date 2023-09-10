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
	void SurfaceReconstruction(int res, Vec3d minPoint, Vec3d maxPoint);  // 曲面重建
	double* resultValues;  // 曲面重建结果((res+1)^3)
private:
	int maxk;  // 最大近邻个数
	double r;//平方支撑半径
	double s;  // CSRBF参数
	double mu;
	int res;  // 曲面重建精度
	Vec3d minPoint, maxPoint;  // 重建范围
	double c;

	int n;  // 点云内点个数
	Vec3d* points;  // 点云坐标数组(生命周期由外部控制)
	Vec3d* normals;  // 点云法向量数组(生命周期由外部控制)
	Octree* octree;  // 点云(生命周期由外部控制)

	Mat3* orthbases;  // 点云各点处局部坐标系的基
	Vec3d* hparas;  // 点云各点对应二次曲面参数
	double* lserror;  // 点云各点二次曲面拟合的最小二乘误差
	double* approxArgs;  // 点云各点对应的拟合系数
	double* cis;  // 点云各点对应的Ci
	double boxlength;


	Nearest* nearest; //近邻信息，大小为n
	double* lambda;
	double* nearestdistance;
	double* mnearestdistance;
	int averagecount;

	inline double RBF(const double& squaredDistanceRatio, const double& approxArg) const;  // 拟合系数控制的的RBF距离比基函数
	inline void CSRBF(const std::vector<int>* Indices, const std::vector<double>* squaredDists, double* values,int k) const;  // 计算给定距离场(由小到大)的CSRBF函数值表
	inline double CalHvalue(const Vec3d& point, const int& i) const;  // 计算给定点在点云第i个点对应的局部二次曲面下的对应h值

	void FindNearestNeighborsAndCalculateHparas(); // 计算点云的近邻信息及每个点的二次曲面参数
	void CalculateApproxArgs();  // 计算各个点的拟合系数
	void CalculateCis();  // 为点云计算每个点的Ci值
	void CalculateReconstructionValues();  // 计算曲面重建结果
	void DeleteNearest();
	void CaculateInditor();

	
};



void IRBF::DeleteNearest()
{
	for(int i=0;i<n;i++)
		nearest[i].Deletenearest();
}

IRBF::IRBF(double r, double s, int n, double mu,Vec3d* points, Vec3d* normals, Octree* octree,double c) {
	// 导入参数
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
	// 初始化曲面重建所需的信息
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
		// 计算该点处局部坐标系的基
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
		// 使用最小二乘法计算该点处二次曲面系数
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
		hparas[i] = tempmat.SolveLinerEquation(tempvec);  // 求解二次曲面系数

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


		// 计算最小二乘法拟合误差
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
	// 求mink和max
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
		double* CSRBFValues = new double[maxk];  // 暂存CSRBF函数值表
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
	int legalCount = 0;  // 识别该点所处的数据点邻域个数
	for (int i = 0; i < findCount; i++) {
		if ((*squaredDists)[i] >= r) {
			continue;  // 该点不在数据点支撑半径内则跳过
		}
		legalCount++;
		t_approx = approxArgs[(*indices)[i]];
		if (t_approx < EPS && (*squaredDists)[i] < EPS) {
			//approxArgs[(*indices)[i]<EPS表示该点为控制点
			chvsum = 0;
			vsum = 1;
			break;
		}
		rbfValue = RBF( lambda[(*indices)[i]]*(*squaredDists)[i] / r, t_approx);
		vsum += rbfValue;
		chvsum += (cis[(*indices)[i]] + CalHvalue(point, (*indices)[i])) * rbfValue;
	}
	if (legalCount < 1) {  // 该点该点所处的数据点邻域个数过少
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
	// 进度条参数
	int progressbar_step = res / 50;
	// 计时
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
				// 计算该点的拟合函数值
				double chvsum = 0, vsum = 0, rbfValue;
				int legalCount = 0;  // 识别该点所处的数据点邻域个数
				for (int i = 0; i < findCount; i++) {
					if ((*squaredDists)[i] >= r) {
						continue;  // 该点不在数据点支撑半径内则跳过
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
				if (legalCount < 1) {  // 该点该点所处的数据点邻域个数过少
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

