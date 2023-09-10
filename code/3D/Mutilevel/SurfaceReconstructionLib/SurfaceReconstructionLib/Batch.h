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
#include"Cube.h"
#include "Surface.h"
#include "rply.h"
#include "NoiseLoad.h"

namespace
{
	bool dotet = false;//使用cube
}



class Batch {
public:
	Batch(const std::string& name,int storageway);
	~Batch();
	bool Execute();
private:
	int autolevel;
	int storageway;
	int level;
	double r;
	double mu;
	Vec3d minPoint, maxPoint;
	double r_size; int bounds;
	double boxlength;
	int res;
	double s;
	std::string name;
	std::string savename;
	FileLoader* loader;
	NoiseLoader* noiseloader;
	ResultSaver* saver;
	Isurface* isurface;
	bool dealnoise;
	Liusurface* liusurface;
	void SaveasIply(std::string savename);
	void SaveasLiuply(std::string savename);
};

Batch::Batch(const std::string& name,int storageway) {
	this->storageway = storageway;
	this->name = name;
	loader = new FileLoader();
	noiseloader = new NoiseLoader();
	saver = new ResultSaver();
	isurface = NULL;
	liusurface = NULL;
	dealnoise = false;
}

Batch::~Batch() {
	if (loader != NULL)
		delete loader;
	if (saver != NULL)
		delete saver;
	if (noiseloader != NULL)
		delete noiseloader;
}

void Batch::SaveasIply(std::string savename)
{
	bounds = 300;
	printf("Reconstructing...\n");
	clock_t start = clock();
	Polygonizer pol(isurface, r_size, bounds, (maxPoint - minPoint).x, (maxPoint - minPoint).y, (maxPoint - minPoint).z); //隐式曲面、网格大小、隐式曲面分量距离限制
	pol.march(dotet, (maxPoint.x + minPoint.x) / 2, (maxPoint.y + minPoint.y) / 2, (maxPoint.z + minPoint.z) / 2);
	double cost = (double)(clock() - start) / 1000;
	printf("Finished in %.3lf seconds.\n", cost);
	p_ply oply;
	int V_N = pol.no_vertices(), F_N = pol.no_triangles();
	oply = ply_create(("results/" + savename + "_irbf.ply").c_str(), PLY_ASCII, NULL);
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
	bounds = 300;
	printf("Reconstructing...\n");
	clock_t start = clock();
	Polygonizer pol(liusurface, r_size, bounds, (maxPoint - minPoint).x, (maxPoint - minPoint).y, (maxPoint - minPoint).z); //隐式曲面、网格大小、隐式曲面分量距离限制
	pol.march(dotet, (maxPoint.x + minPoint.x) / 2, (maxPoint.y + minPoint.y) / 2, (maxPoint.z + minPoint.z) / 2);
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
			r = 0.75 * boxlength;
			autolevel = ceil(-log2(this->loader->rhohat / (2 * r)));
			std::cout << "autolevel:" << autolevel << std::endl;
			//std::cout << "auto:" << 2* this->loader->rhohat/r << std::endl;
		}
		else if (para == "deal_noise")
		{
			file >> para;
			noiseloader->Load(para);
			dealnoise = true;

		}
		else if (para == "IRBF") {
			file >> para;
			if (para == "Init") {
				file >> level;
				level = autolevel;
				 double mu;
				file >> s; file >> res;
				mu = 0.1 * r * r;
				printf("IRBF Init: r=%f, s=%f\n", r, s);
				isurface = new Isurface(level,r, s,mu,res,loader->cubetree->root,minPoint,maxPoint,loader->points,loader->normals,loader->n,
					noiseloader->noise,noiseloader->n,dealnoise);
				/*double maxe = 0; double meane = 0;
				for (int i = 0; i < loader->n; i++)
				{
					double temp = abs(isurface->eval(loader->points[i].x, loader->points[i].y, loader->points[i].z));
					maxe = maxe > temp ? maxe : temp;
					meane += temp;
				}
				meane = meane / loader->n;
				std::cout << "最大误差：" << maxe << std::endl;
				std::cout << "平均误差：" << meane << std::endl;*/

			}
			else if (para == "Reconstruct") {
				//Vec3d resminPoint, resmaxPoint;
				file >> r_size;
				//std::cout << iRBF->eval(-0.063712, 0.044575, 0.033687);
				//file >> resminPoint.x >> resminPoint.y >> resminPoint.z;
				//file >> resmaxPoint.x >> resmaxPoint.y >> resmaxPoint.z;
				file >> savename;
				if (storageway == 1)
				{
					savename = savename + "_level" + std::to_string(level) + "_s" + std::to_string(s) + "_res" + std::to_string(res) + "_irbf";
					std::cout << "Saving result as results/" << para << ".mat...\n";
					isurface->Reconstruction();
					saver->Save(res, isurface->resultvalues, savename);
					printf("Save finished.\n");
				}
				else {
					savename = savename + "_level" + std::to_string(level) + "_s" + std::to_string(s) + "_size" + std::to_string(r_size) + "_irbf";
					SaveasIply(savename);
				}
				
			}
			else if (para == "Release") {
				delete isurface;
				isurface = NULL;
				printf("Released IRBF resources.\n");
			}
		}

		else if (para == "LiuRBF") {
			file >> para;
			if (para == "Init") {
				file >> level; file >> res;
				level = autolevel;
				double mu;
				mu = 0.1 * r * r;
				printf("LiuRBF Init: r=%f\n", r);
				liusurface = new Liusurface(level, r, mu,res, loader->cubetree->root, minPoint, maxPoint);
				double maxe = 0; double meane = 0;
				/*for (int i = 0; i < loader->n; i++)
				{
					double temp = abs(liusurface->eval(loader->points[i].x, loader->points[i].y, loader->points[i].z));
					maxe = maxe > temp ? maxe : temp;
					meane += temp;
				}
				meane = meane / loader->n;
				std::cout << "最大误差：" << maxe << std::endl;
				std::cout << "平均误差：" << meane << std::endl;*/
			}
			else if (para == "Reconstruct") {
				file >> r_size;
				//std::cout << iRBF->eval(-0.063712, 0.044575, 0.033687);
				file >> savename;
				if (storageway == 1)
				{
					savename = savename + "_level" + std::to_string(level) + "_s" + std::to_string(s) + "_res" + std::to_string(res)+ "_liurbf";
					std::cout << "Saving result as results/" << para << ".mat...\n";
					//liusurface->Reconstruction();
					//saver->Save(res, liusurface->resultvalues, savename);
					printf("Save finished.\n");
				}
				else {
					savename = savename + "_level" + std::to_string(level) + "_s" + std::to_string(s) + "_size" + std::to_string(r_size) + "_liurbf";
					SaveasLiuply(savename);
				}
			}
			else if (para == "Release") {
				delete liusurface;
				liusurface = NULL;
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