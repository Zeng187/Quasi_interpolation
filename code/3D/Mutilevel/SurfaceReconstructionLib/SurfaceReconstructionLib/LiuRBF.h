#pragma once
#include "Vec3d.h"
#include "Mat3.h"
#include "Octree.h"
#include <vector>
#include <cmath>
#include <time.h>
#include <float.h>
#include"Nearest.h"

#define EPS 1e-10


class LiuRBF {
public:
	LiuRBF(double r,double mu, int n,int res,double T, Vec3d* points, Vec3d* normals, Octree* octree);
	~LiuRBF();
	double* resultValues;  // 曲面重建结果((res+1)^3)
	void SurfaceReconstruction(int res, Vec3d minPoint, Vec3d maxPoint);  // 曲面重建
	void CalculateReconstructionValues(Vec3d* Dpoints, int tn, double* gvalue);
	void RenewCis(double* gvalue, int tn);  // 在层次网格中，使第k层的每个点的Ci值加上上一层的计算结果
	inline double eval(double x, double y, double z);
private:
	int maxk;  // 最大近邻个数
	double r;//平方支撑半径
	int res;  // 曲面重建精度
	Vec3d minPoint, maxPoint;  // 重建范围
	double mu;//形状参数
	int n;  // 点云内点个数
	Vec3d* points;  // 点云坐标数组(生命周期由外部控制)
	Vec3d* normals;  // 点云法向量数组(生命周期由外部控制)
	Octree* octree;  // 点云(生命周期由外部控制)

	Mat3* orthbases;  // 点云各点处局部坐标系的基
	Vec3d* hparas;  // 点云各点对应二次曲面参数
	double* cis;  // 点云各点对应的Ci
	double* lambda;
	double T;

	Nearest* nearest; //近邻信息，大小为n
	double* nearestdistance;
	int adaptive_count;
	double maxsquareddistance;

	inline double RBF(const double& squaredDistanceRatio) const;  // 拟合系数控制的的RBF距离比基函数
	inline void CSRBF(const std::vector<int>* Indices, const std::vector<double>* squaredDists, double* values, int k) const;  // 计算给定距离场(由小到大)的CSRBF函数值表
	inline double CalHvalue(const Vec3d& point, const int& i) const;  // 计算给定点在点云第i个点对应的局部二次曲面下的对应h值

	void FindNearestNeighborsAndCalculateHparas();  // 计算点云的k近邻信息及每个点的二次曲面参数
	void CalculateCis();  // 为点云计算每个点的Ci值
	void CalculateReconstructionValues();  // 计算曲面重建结果
	void DeleteNearest();
};

LiuRBF::LiuRBF(double r,double mu, int n, int res,double T,Vec3d* points, Vec3d* normals, Octree* octree) {
	// 导入参数
	this->r = r;
	this->n = n;
	this->points = points;
	this->normals = normals;
	this->octree = octree;
	this->mu = mu;
	this->res = res;
	this->T = T;
	orthbases = new Mat3[n];
	hparas = new Vec3d[n];
	cis = new double[n];
	nearest = new Nearest[n];
	nearestdistance = new double[n];
	lambda = new double[n];
	resultValues = NULL;
	adaptive_count =16 < n ? 16 : n;
	//adaptive_count =0;
	maxsquareddistance = 0;
	// 初始化曲面重建所需的信息
	//printf("Precalculating info for surface reconstruction...\n");
	//clock_t start = clock();

	FindNearestNeighborsAndCalculateHparas();
	CalculateCis();


	//double cost = (double)(clock() - start) / 1000;
	//printf("Finished in %.3lf seconds.\n", cost);
}

LiuRBF::~LiuRBF() {
	delete[] orthbases;
	delete[] hparas;
	delete[] cis;
	delete[] nearestdistance;
	delete[] lambda;
	DeleteNearest();
	if (resultValues != NULL) { delete[] resultValues; }
}


inline double LiuRBF::RBF(const double& r) const {
	double sqrtr = sqrt(r); double temp = (1 - sqrtr) * (1 - sqrtr);
	return temp*temp * (4 * sqrtr + 1);
}

inline void LiuRBF::CSRBF(const std::vector<int>* Indices, const std::vector<double>* squaredDists, double* values, int k) const {
	double sumVal = 0;
	for (int i = 0; i < k; i++) {
		values[i] = RBF(lambda[(*Indices)[i]] * (*squaredDists)[i] /  nearest[(*Indices)[i]].r);
		sumVal += values[i];
	}
	for (int i = 0; i < k; i++) {
		values[i] /= sumVal;
	}
}

void LiuRBF::SurfaceReconstruction(int res, Vec3d minPoint, Vec3d maxPoint) {
	this->res = res;
	this->minPoint = minPoint;
	this->maxPoint = maxPoint;
	if (resultValues != NULL) { delete[] resultValues; }
	resultValues = new double[(res + 1) * (res + 1) * (res + 1)];
	CalculateReconstructionValues();
}

inline double LiuRBF::CalHvalue(const Vec3d& point, const int& i) const {
	Vec3d transformedLoc = orthbases[i].BaseTransform(point - points[i]);
	double surfaceW =
		hparas[i].x * transformedLoc.x * transformedLoc.x +
		hparas[i].y * transformedLoc.x * transformedLoc.y +
		hparas[i].z * transformedLoc.y * transformedLoc.y;
	return (transformedLoc.z - surfaceW);
}

void LiuRBF::FindNearestNeighborsAndCalculateHparas() {
	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		Nearest t_nearest = nearest[i];
		t_nearest.count = octree->FindNearestNeighborsBySquaredRadius(points[i], r, t_nearest.indices, t_nearest.squaredDists);
		double r1 = r;
		while (t_nearest.count < adaptive_count)
		{
			r1 = 1.1 * r1;
			t_nearest.count = octree->FindNearestNeighborsBySquaredRadius(points[i], r1, t_nearest.indices, t_nearest.squaredDists);
		}
		nearestdistance[i] = (*t_nearest.squaredDists)[1];
		nearest[i].count = t_nearest.count;
		nearest[i].r = r1;
		//lambda[i] = sqrt(mu / nearestdistance[i]);
		lambda[i] = 1;//不设形状参数
		maxk = t_nearest.count > maxk ? t_nearest.count : maxk;
		maxsquareddistance = maxsquareddistance > r1 ? maxsquareddistance : r1;
		// 计算该点处局部坐标系的基
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
			tempmat.u.x += transformedLoc.x * transformedLoc.x * transformedLoc.x * transformedLoc.x;
			tempmat.v.x += transformedLoc.x * transformedLoc.x * transformedLoc.x * transformedLoc.y;
			tempmat.w.x += transformedLoc.x * transformedLoc.x * transformedLoc.y * transformedLoc.y;
			tempmat.w.y += transformedLoc.x * transformedLoc.y * transformedLoc.y * transformedLoc.y;
			tempmat.w.z += transformedLoc.y * transformedLoc.y * transformedLoc.y * transformedLoc.y;
			tempvec.x += transformedLoc.x * transformedLoc.x * transformedLoc.z;
			tempvec.y += transformedLoc.x * transformedLoc.y * transformedLoc.z;
			tempvec.z += transformedLoc.y * transformedLoc.y * transformedLoc.z;
		}
		tempmat.u.y = tempmat.v.x;
		tempmat.u.z = tempmat.v.y = tempmat.w.x;
		tempmat.v.z = tempmat.w.y;
		hparas[i] = tempmat.SolveLinerEquation(tempvec);  // 求解二次曲面系数
	}

}

void LiuRBF::CalculateCis() {
	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		double* CSRBFValues = new double[maxk];  // 暂存CSRBF函数值表
		Nearest t_nearest = nearest[i];
		CSRBF(t_nearest.indices,t_nearest.squaredDists, CSRBFValues, t_nearest.count);
		double c = 0;
		for (int j = 0; j < t_nearest.count; j++) {
			c -= CalHvalue(points[i], (*t_nearest.indices)[j]) * CSRBFValues[j];
		}
		cis[i] = c/(1+T);
		//cis[i] = c;
		delete[] CSRBFValues;
	}
	
}


void LiuRBF::CalculateReconstructionValues() {
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
		std::vector<int>* indices = new std::vector<int>;
		std::vector<double>* squaredDists = new std::vector<double>;
		for (int ix = 0; ix <= res; ix++) {
			for (int iy = 0; iy <= res; iy++) {
				Vec3d point(ix, iy, iz);
				point.ApplyMultiplication(stepVec);
				point += minPoint;
				int findCount = octree->FindNearestNeighborsBySquaredRadius(point, maxsquareddistance, indices, squaredDists);
				int index = iz * sqcount + ix * count + iy;
				// 计算该点的拟合函数值
				double chvsum = 0, vsum = 0, rbfValue;
				int legalCount = 0;  // 识别该点所处的数据点邻域个数
				for (int i = 0; i < findCount; i++) {
					if ((*squaredDists)[i] >= nearest[(*indices)[i]].r) {
						continue;  // 该点不在数据点支撑半径内则跳过
					}
					legalCount++;
					rbfValue = RBF(lambda[(*indices)[i]] * (*squaredDists)[i] / nearest[(*indices)[i]].r);
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

void LiuRBF::CalculateReconstructionValues(Vec3d* Dpoints, int tn, double* gvalue)
{

	Vec3d point; int legalCount; double chvsum, vsum, rbfValue; int findCount;
	std::vector<int>* indices = new std::vector<int>;
	std::vector<double>* squaredDists = new std::vector<double>;
	for (int i = 0; i < tn; i++)
	{
		point = Dpoints[i];
		chvsum = 0, vsum = 0;
		legalCount = 0;  // 识别该点所处的数据点邻域个数
		findCount = octree->FindNearestNeighborsBySquaredRadius(Dpoints[i], maxsquareddistance, indices, squaredDists);
		for (int j = 0; j < findCount; j++) {
			if ((*squaredDists)[j] >= nearest[(*indices)[j]].r) {
				continue;  // 该点不在数据点支撑半径内则跳过
			}
			legalCount++;
			rbfValue = RBF(lambda[(*indices)[j]] * (*squaredDists)[j] / nearest[(*indices)[j]].r);
			vsum += rbfValue;
			chvsum += (cis[(*indices)[j]] + CalHvalue(point, (*indices)[j])) * rbfValue;
		}
		if (legalCount < 1) {  // 该点该点所处的数据点邻域个数过少
			gvalue[i] = 32768;
		}
		else {
			gvalue[i] = chvsum / vsum;
		}
	}
	delete indices;
	delete squaredDists;
}

void LiuRBF::RenewCis(double* gvalue, int tn)
{
	for (int i = 0; i < tn; i++)
	{
		this->cis[i] = -gvalue[i] + this->cis[i];
		this->cis[i] = this->cis[i] / (1 + T);
	}
}


void LiuRBF::DeleteNearest()
{
	for (int i = 0; i < n; i++)
		nearest[i].Deletenearest();
}

double LiuRBF::eval(double x, double y, double z)
{
	Vec3d point = Vec3d(x, y, z);
	double chvsum = 0, vsum = 0, rbfValue;
	int legalCount = 0;  // 识别该点所处的数据点邻域个数
	std::vector<int>* indices = new std::vector<int>;
	std::vector<double>* squaredDists = new std::vector<double>;
	int findCount = octree->FindNearestNeighborsBySquaredRadius(point, maxsquareddistance, indices, squaredDists);
	for (int j = 0; j < findCount; j++) {
		if ((*squaredDists)[j] >= nearest[(*indices)[j]].r) {
			continue;  // 该点不在数据点支撑半径内则跳过
		}
		legalCount++;
		rbfValue = RBF(lambda[(*indices)[j]] * (*squaredDists)[j] / nearest[(*indices)[j]].r);
		vsum += rbfValue;
		chvsum += (cis[(*indices)[j]] + CalHvalue(point, (*indices)[j])) * rbfValue;
	}
	if (legalCount < 1) {  // 该点该点所处的数据点邻域个数过少
		return(32768);
	}
	else {
		return(chvsum / vsum);
	}
	delete indices;
	delete squaredDists;
}