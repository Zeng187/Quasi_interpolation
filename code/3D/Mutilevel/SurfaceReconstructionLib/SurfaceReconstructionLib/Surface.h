#pragma once
#include"polygonizer.h"
#include"IRBF.h"
#include"LiuRBF.h"
#include"Cube.h"
#include "Octree.h"
#include <time.h>

class Isurface :public ImplicitFunction
{
public:
	int level;
	double* resultvalues;
	Isurface(int level, double r, double s, double mu, int res, CubeNode* root, Vec3d minPoint, Vec3d maxPoint, Vec3d* points, Vec3d* normals, int n,
		int *noise,int noisecount,bool dealnoise);
	inline virtual double eval(double x, double y, double z);
	~Isurface();
	void Reconstruction();
private:
	double r;
	double s;
	double mu;
	int res;
	Vec3d* points, * normals; int n;
	IRBF** iRBF;
	CubeNode* h;
	int* tcount;//记录每一层的数量
	Vec3d** Dpoints;//每一层的点
	Vec3d** Dnormals;//每一层点对应的法向量
	Octree** octree;//每一层由点云产生的八叉树结构
	int* noise;
	int noisecount;
	Vec3d minPoint;
	Vec3d maxPoint;
	double T;
	bool dealnoise;


};

Isurface::Isurface(int level, double r, double s, double mu, int res, CubeNode* root, Vec3d minPoint, Vec3d maxPoint,  Vec3d* points, Vec3d* normals, int n,
	int* noise, int noisecount,bool dealnoise)
{
	this->level = level;
	this->r = r;
	this->s = s;
	this->mu = mu;
	this->res = res;
	this->points = points; this->normals = normals; this->n = n;
	h = new CubeNode();
	tcount = new int[level];
	iRBF = new IRBF * [level];
	Dpoints = new Vec3d * [level];
	Dnormals = new Vec3d * [level];
	octree = new Octree * [level];
	h->next = root;
	this->minPoint = minPoint;
	this->maxPoint = maxPoint;
	this->resultvalues = new double[(res + 1) * (res + 1) * (res + 1)];
	this->noise = noise;
	this->noisecount = noisecount;
	this->dealnoise = dealnoise;

	double** gvalue = new double* [level];
	double** t_gvalue = new double* [level];


	clock_t start = clock();
	int i, j, k;
	double temp_r = r; double temp_mu = mu;
	double* c = new double[level];

	for (i = 0; i < level; i++)
	{
		c[i] = 1.0/((i+1)*(i+1));
	}
	//构造
	double T; 
	//temp_mu = 0.01;
	for (i = 0; i < level; i++)
	{
		T = 0.25*(i+1);
		std::cout << temp_mu << std::endl;
		if (false)
		{
			tcount[i] = n;
			Dpoints[i] = points;
			Dnormals[i] = normals;
			octree[i] = new Octree(minPoint, maxPoint, 8);
			octree[i]->AddPoints(points, n);
			iRBF[i] = new IRBF(temp_r * temp_r, s, temp_mu, tcount[i], res, c[i], T, Dpoints[i], Dnormals[i], octree[i], true, noise, noisecount, dealnoise);

			if (i != 0)
			{

				//std::cout << i << std::endl;
				gvalue[i] = new double[tcount[i]];
				//初始化gvalue的值
#pragma omp parallel for
				for (j = 0; j < tcount[i]; j++)
				{
					gvalue[i][j] = 0;
				}
#pragma omp parallel for
				for (j = 0; j < i; j++)
				{
					t_gvalue[j] = new double[tcount[i]];
					iRBF[j]->CalculateReconstructionValues(Dpoints[i], tcount[i], t_gvalue[j]);
					for (k = 0; k < tcount[i]; k++)
					{
						gvalue[i][k] += t_gvalue[j][k];
					}
					delete[] t_gvalue[j];
				}

				iRBF[i]->RenewCis(gvalue[i], tcount[i]);

			}

		}
		else
		{
			tcount[i] = h->ConstructCubeNodelink(h);
			Dpoints[i] = new Vec3d[tcount[i]];
			Dnormals[i] = new Vec3d[tcount[i]];
			h->Collect(Dpoints[i], Dnormals[i], h, tcount[i]);
			octree[i] = new Octree(minPoint, maxPoint, 8);
			octree[i]->AddPoints(Dpoints[i], tcount[i]);
			iRBF[i] = new IRBF(temp_r * temp_r, s, temp_mu, tcount[i], res, c[i],T, Dpoints[i], Dnormals[i], octree[i], false, noise, noisecount, false);

			if (i != 0)
			{

				//std::cout << i << std::endl;
				gvalue[i] = new double[tcount[i]];
				//初始化gvalue的值
#pragma omp parallel for
				for (j = 0; j < tcount[i]; j++)
				{
					gvalue[i][j] = 0;
				}
#pragma omp parallel for
				for (j = 0; j < i; j++)
				{
					t_gvalue[j] = new double[tcount[i]];
					iRBF[j]->CalculateReconstructionValues(Dpoints[i], tcount[i], t_gvalue[j]);
					for (k = 0; k < tcount[i]; k++)
					{
						gvalue[i][k] += t_gvalue[j][k];
					}
					delete[] t_gvalue[j];
				}


				iRBF[i]->RenewCis(gvalue[i], tcount[i]);
			}

		}

		temp_r = temp_r / 2;
		temp_mu = temp_mu / 4;
	}
	iRBF[0]->writeinformation();
	double cost = (double)(clock() - start) / 1000;
	printf("Finished in %.3lf seconds.\n", cost);
	if (level != 1)
	{
		for (j = 1; j < level; j++)
		{
			delete[] gvalue[j];
		}
	}
	delete[] gvalue;
	delete[] t_gvalue;
	delete[] c;
}

Isurface::~Isurface()
{
	for (int i = 0; i < level; i++)
	{
		delete iRBF[i];
		delete octree[i];
		if (i < level - 1)
		{
			delete[] Dpoints[i];
			delete[] Dnormals[i];
		}

	}
	delete[] Dpoints;
	delete[] Dnormals;
	if (resultvalues != NULL)
		delete[] resultvalues;
}


double Isurface::eval(double x, double y, double z)
{
	double value = 0, temp;
	bool inside = false;
	for (int i = 0; i < level; i++)
	{
		temp = iRBF[i]->eval(x, y, z);
		if (temp < 32768)
		{
			inside = true;
			value += temp;
		}

	}
	if (!inside)
		value = 32768;
	return value;
}



void Isurface::Reconstruction()
{
	double** values = new double* [level]; int i, j;
	bool* inside = new bool[(res + 1) * (res + 1) * (res + 1)];
	for (i = 0; i < level; i++)
	{
		iRBF[i]->SurfaceReconstruction(res, minPoint, maxPoint);
		values[i] = iRBF[i]->resultValues;
		iRBF[i]->resultValues = NULL;
	}
	for (j = 0; j < (res + 1) * (res + 1) * (res + 1); j++)
	{
		inside[j] = false;
		for (i = 0; i < level; i++)
		{
			if (values[i][j] < 32768)
			{
				//正常在范围内的
				resultvalues[j] += values[i][j];
				inside[j] = true;
			}
			//否则，不在这一层的任何一个函数支撑范围内

		}
		if (!inside[j])
		{
			//如果点不在任何一层的任何一个函数支撑范围内，则赋予其充分大的值
			resultvalues[j] = 32768;
		}
	}
	for (j = 0; j < level; j++)
	{
		delete[] values[j];
	}
	delete[] values;
	delete[] inside;

}









class Liusurface :public ImplicitFunction
{
public:
	int level;
	double* resultvalues;
	Liusurface(int level, double r, double mu, int res, CubeNode* root, Vec3d minPoint, Vec3d maxPoint);
	inline virtual double eval(double x, double y, double z);
	~Liusurface();
	void Reconstruction();
private:
	double r;
	double s;
	double mu;
	int res;
	LiuRBF** liuRBF;
	CubeNode* h;
	int* tcount;//记录每一层的数量
	Vec3d** Dpoints;//每一层的点
	Vec3d** Dnormals;//每一层点对应的法向量
	Octree** octree;//每一层由点云产生的八叉树结构
	Vec3d minPoint;
	Vec3d maxPoint;


};




Liusurface::Liusurface(int level, double r, double mu, int res, CubeNode* root, Vec3d minPoint, Vec3d maxPoint)
{
	this->level = level;
	this->r = r;
	this->mu = mu;
	this->res = res;
	liuRBF = new LiuRBF * [level];
	h = new CubeNode();
	tcount = new int[level];
	Dpoints = new Vec3d * [level];
	Dnormals = new Vec3d * [level];
	octree = new Octree * [level];
	h->next = root;
	this->minPoint = minPoint;
	this->maxPoint = maxPoint;
	this->resultvalues = new double[(res + 1) * (res + 1) * (res + 1)];

	double** gvalue = new double* [level];
	double** t_gvalue = new double* [level];



	clock_t start = clock();
	int i, j, k;
	double temp_r = r; double temp_mu = mu;
	//构造
	double T;
	for (i = 0; i < level; i++)
	{
		tcount[i] = h->ConstructCubeNodelink(h);
		Dpoints[i] = new Vec3d[tcount[i]];
		Dnormals[i] = new Vec3d[tcount[i]];
		h->Collect(Dpoints[i], Dnormals[i], h, tcount[i]);
		octree[i] = new Octree(minPoint, maxPoint, 8);
		octree[i]->AddPoints(Dpoints[i], tcount[i]);
		T = 0.25 * (1 + i);
		liuRBF[i] = new LiuRBF(temp_r * temp_r, temp_mu, tcount[i], res, T,Dpoints[i], Dnormals[i], octree[i]);
		std::cout << i << std::endl;
		if (i != 0)
		{
			gvalue[i] = new double[tcount[i]];
			//初始化gvalue的值
			for (j = 0; j < tcount[i]; j++)
			{
				gvalue[i][j] = 0;
			}
#pragma omp parallel for
			for (j = 0; j < tcount[i]; j++)
			{
				gvalue[i][j] = 0;
			}
#pragma omp parallel for
			for (j = 0; j < i; j++)
			{
				t_gvalue[j] = new double[tcount[i]];
				liuRBF[j]->CalculateReconstructionValues(Dpoints[i], tcount[i], t_gvalue[j]);
				for (k = 0; k < tcount[i]; k++)
				{
					gvalue[i][k] += t_gvalue[j][k];
				}
				delete[] t_gvalue[j];
			}

			liuRBF[i]->RenewCis(gvalue[i], tcount[i]);
		}

		temp_r = temp_r / 2;
		temp_mu = temp_mu / 4;
	}
	double cost = (double)(clock() - start) / 1000;
	printf("Finished in %.3lf seconds.\n", cost);
	if (level != 1)
	{
		for (j = 1; j < level; j++)
		{
			delete[] gvalue[j];
		}
	}
	delete[] gvalue;
	delete[] t_gvalue;
}


Liusurface::~Liusurface()
{
	for (int i = 0; i < level; i++)
	{
		delete liuRBF[i];
		delete octree[i];
		delete[] Dpoints[i];
		delete[] Dnormals[i];
	}
	delete[] Dpoints;
	delete[] Dnormals;
	if (resultvalues != NULL)
		delete[] resultvalues;
}

double Liusurface::eval(double x, double y, double z)
{
	double value = 0, temp;
	bool inside = false;
	for (int i = 0; i < level; i++)
	{
		temp = liuRBF[i]->eval(x, y, z);
		if (temp < 32768)
		{
			inside = true;
			value += temp;
		}

	}
	if (!inside)
		value = 32768;
	return value;
}


void Liusurface::Reconstruction()
{
	double** values = new double* [level]; int i, j;
	bool* inside = new bool[(res + 1) * (res + 1) * (res + 1)];
	for (i = 0; i < level; i++)
	{
		liuRBF[i]->SurfaceReconstruction(res, minPoint, maxPoint);
		values[i] = liuRBF[i]->resultValues;
		liuRBF[i]->resultValues = NULL;
	}
	for (j = 0; j < (res + 1) * (res + 1) * (res + 1); j++)
	{
		inside[j] = false;
		for (i = 0; i < level; i++)
		{
			if (values[i][j] < 32768)
			{
				//正常在范围内的
				resultvalues[j] += values[i][j];
				inside[j] = true;
			}
			//否则，不在这一层的任何一个函数支撑范围内

		}
		if (!inside[j])
		{
			//如果点不在任何一层的任何一个函数支撑范围内，则赋予其充分大的值
			resultvalues[j] = 32768;
		}
	}
	for (j = 0; j < level; j++)
	{
		delete[] values[j];
	}
	delete[] values;
	delete[] inside;

}