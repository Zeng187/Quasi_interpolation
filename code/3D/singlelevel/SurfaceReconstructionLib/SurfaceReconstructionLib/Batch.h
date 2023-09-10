#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <cstdio>
#include "Octree.h"
#include "Vec3d.h"
#include "Mat3.h"
#include "IRBF.h"
#include "LiuRBF.h"
#include "FileLoader.h"
#include "ResultSaver.h"
#include <cstdio>
#include <string>
#include <stdio.h> 
#include "polygonizer.h"
#include "rply.h"
#include <ctime>


namespace
{
	bool dotet = true;//使用cube
}


class Batch {
public:
	Batch(const std::string& name, int storageway);
	void SaveasIply(std::string savename);
	void SaveasLiuply(std::string savename);
	~Batch();
	bool Execute();
private:
	double r;
	double autor;
	double s; int res;
	double c;
	int storageway;
	Vec3d minPoint, maxPoint;
	double r_size; int bounds;
	double boxlength;
	std::string name;
	std::string savename;
	FileLoader* loader;
	ResultSaver* saver;
	IRBF* iRBF;
	LiuRBF* liuRBF;
};

Batch::Batch(const std::string& name, int storageway) {
	this->storageway = storageway;
	this->name = name;
	loader = new FileLoader();
	saver = new ResultSaver();
}


Batch::~Batch() {
	delete loader;
	delete saver;
	if (iRBF) { delete iRBF; }
	if (liuRBF) { delete liuRBF; }
}



void Batch::SaveasIply(std::string savename)
{
	bounds = 2 / r_size > 200 ? 2 / r_size : 200;
	printf("Reconstructing...\n");
	clock_t start = clock();
	Polygonizer pol(iRBF, r_size, bounds, (maxPoint - minPoint).x, (maxPoint - minPoint).y, (maxPoint - minPoint).z); //隐式曲面、网格大小、隐式曲面分量距离限制
	pol.march(dotet, (maxPoint.x + minPoint.x) / 2, (maxPoint.y + minPoint.y) / 2, (maxPoint.z + minPoint.z) / 2);
	double cost = (double)(clock() - start) / 1000;
	printf("Finished in %.3lf seconds.\n", cost);
	p_ply oply;
	int V_N = pol.no_vertices(), F_N = pol.no_triangles();
	oply = ply_create(("results/"+savename+"_irbf.ply").c_str(), PLY_ASCII, NULL);
	//添加vertex元素
	ply_add_element(oply, "vertex", V_N);
	//添加vertex的属性
	ply_add_property(oply, "x", PLY_DOUBLE, PLY_LIST, PLY_LIST);
	ply_add_property(oply, "y", PLY_DOUBLE, PLY_LIST, PLY_LIST);
	ply_add_property(oply, "z", PLY_DOUBLE, PLY_LIST, PLY_LIST);
	//添加face元素
	ply_add_element(oply, "face", F_N);
	//添加face属性
	ply_add_property(oply, "vertex_indices", PLY_LIST, PLY_INT, PLY_INT);
	//将header部分写入ply文件
	ply_write_header(oply);
	//将数据项写入ply文件
	TRIANGLE t; int temp = 0;
	for (int i = 0; i < V_N; i++)
	{
		ply_write(oply, (pol.gvertices)[i].x);
		ply_write(oply, (pol.gvertices)[i].y);
		ply_write(oply, (pol.gvertices)[i].z);
	}
	for (int j = 0; j < F_N; j++)
	{
		ply_write(oply, 3);
		ply_write(oply, (pol.gtriangles)[j].v0);
		ply_write(oply, (pol.gtriangles)[j].v1);
		ply_write(oply, (pol.gtriangles)[j].v2);
	}
	ply_close(oply);
	printf("Saveing finished!\n");
}


void Batch::SaveasLiuply(std::string savename)
{
	bounds = 2 / r_size > 200 ? 2 / r_size : 200;
	printf("Reconstructing...\n");
	clock_t start = clock();
	Polygonizer pol(liuRBF, r_size, bounds, (maxPoint - minPoint).x, (maxPoint - minPoint).y, (maxPoint - minPoint).z); //隐式曲面、网格大小、隐式曲面分量距离限制
	pol.march(dotet, (maxPoint.x + minPoint.x) / 2, (maxPoint.y + minPoint.y) / 2, (maxPoint.z + minPoint.z) /  2);
	double cost = (double)(clock() - start) / 1000;
	printf("Finished in %.3lf seconds.\n", cost);
	p_ply oply;
	int V_N = pol.no_vertices(), F_N = pol.no_triangles();
	oply = ply_create(("results/" + savename + "_liurbf.ply").c_str(), PLY_ASCII, NULL);
	//添加vertex元素
	ply_add_element(oply, "vertex", V_N);
	//添加vertex的属性
	ply_add_property(oply, "x", PLY_DOUBLE, PLY_LIST, PLY_LIST);
	ply_add_property(oply, "y", PLY_DOUBLE, PLY_LIST, PLY_LIST);
	ply_add_property(oply, "z", PLY_DOUBLE, PLY_LIST, PLY_LIST);
	//添加face元素
	ply_add_element(oply, "face", F_N);
	//添加face属性
	ply_add_property(oply, "vertex_indices", PLY_LIST, PLY_INT, PLY_INT);
	//将header部分写入ply文件
	ply_write_header(oply);
	//将数据项写入ply文件
	TRIANGLE t; int temp = 0;
	for (int i = 0; i < V_N; i++)
	{
		ply_write(oply, (pol.gvertices)[i].x);
		ply_write(oply, (pol.gvertices)[i].y);
		ply_write(oply, (pol.gvertices)[i].z);
	}
	for (int j = 0; j < F_N; j++)
	{
		ply_write(oply, 3);
		ply_write(oply, (pol.gtriangles)[j].v0);
		ply_write(oply, (pol.gtriangles)[j].v1);
		ply_write(oply, (pol.gtriangles)[j].v2);
	}
	ply_close(oply);
	printf("Save finished!\n");
}



bool Batch::Execute() {
	std::ifstream file;
	file.open(name.c_str());
	if (!file) {
		return false;
	}
	std::string para;
	while (true) {
		file >> para;
		if (para == "Begin") { break; }
	}
	while (true) {
		file >> para;
		if (para == "Load") {
			file >> para;
			loader->Load(para);
			minPoint = loader->minPoint; maxPoint = loader->maxPoint;
			boxlength = (maxPoint - minPoint).Length();
			int leafcount = 0;
			autor = this->loader->octree->Averageleaflength(this->loader->octree->root, leafcount);
			c= 2 * autor / leafcount / boxlength ;
			autor = 0.75*autor / leafcount;

		}

		else if (para == "IRBF") {
			file >> para;
			if (para == "Init") {
				file >> r;
				r = autor;
				double mu = 0.1 * r * r;
				//mu = 0.01;
				//std::cout << "mu:"<<mu << std::endl;
				file >> s; file >> res;
				printf("IRBF Init: r=%f, s=%f\n", r, s);
				iRBF = new IRBF(r*r, s, loader->n,mu, loader->points, loader->normals, loader->octree,c);
				double maxe = 0;double meane = 0;
				for (int i = 0; i < loader->n; i++)
				{
					double temp = abs(iRBF ->eval(loader->points[i].x, loader->points[i].y, loader->points[i].z));
					maxe = maxe > temp ? maxe : temp;
					meane = meane + temp;
				}
				meane = meane / loader->n;
				std::cout << "最大误差：" << maxe << std::endl;
				std::cout << "平均误差：" << meane << std::endl;
			}
			else if (para == "Reconstruct") {
				file >> r_size;
				//std::cout << iRBF->eval(-0.063712, 0.044575, 0.033687);
				file >> savename;
				if (storageway == 1)
				{
					savename = savename  + "_s" + std::to_string(s) + "_res" + std::to_string(res) + "_irbf";
					std::cout << "Saving result as results/" << para << ".mat...\n";
					iRBF->SurfaceReconstruction(res,this->minPoint,this->maxPoint);
					saver->Save(res, iRBF->resultValues, savename);
					printf("Save finished.\n");
				}
				else {
					savename = savename+ "_s" + std::to_string(s) + "_size" + std::to_string(r_size) + "_irbf";
					SaveasIply(savename);
				}
			}
			else if (para == "Release") {
				delete iRBF;
				iRBF = NULL;
				printf("Released IRBF resources.\n");
			}
		}

		else if (para == "LiuRBF") {
			file >> para;
			if (para == "Init") {
				file >> r; file >> res;
				r = autor;
				double mu = 0.1 * r * r;
				//mu = 0.1;
				//std::cout << "mu:" << mu << std::endl;
				printf("LiuRBF Init: r=%f\n", r);
				liuRBF = new LiuRBF(r*r, loader->n,mu, loader->points, loader->normals, loader->octree);
				double maxe = 0; double meane = 0;
				for (int i = 0; i < loader->n; i++)
				{
					double temp = abs(liuRBF->eval(loader->points[i].x, loader->points[i].y, loader->points[i].z));
					maxe = maxe > temp ? maxe : temp;
					meane = meane + temp;
				}
				meane = meane / loader->n;
				std::cout << "最大误差：" << maxe << std::endl;
				std::cout << "平均误差：" << meane << std::endl;
			}
			else if (para == "Reconstruct") {
				file >> r_size;
				file >> savename;
				if (storageway == 1)
				{
					savename = savename  + "_res" + std::to_string(res) + "_liurbf";
					std::cout << "Saving result as results/" << para << ".mat...\n";
					liuRBF->SurfaceReconstruction(res, this->minPoint, this->maxPoint);
					saver->Save(res, liuRBF->resultValues, savename);
					printf("Save finished.\n");
				}
				else {
					savename = savename + "_size" + std::to_string(r_size) + "_liurbf";
					SaveasLiuply(savename);
				}
			}
			else if (para == "Release") {
				delete liuRBF;
				liuRBF = NULL;
				printf("Released LiuRBF resources.\n");
			}
		}

		else if (para == "End") {
			break;
		}
	}
	file.close();
	return true;
}