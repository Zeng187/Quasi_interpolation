#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <cstdio>
#include "Vec3d.h"
#include "Octree.h"

class FileLoader {
public:
	FileLoader();
	~FileLoader();
	bool Load(const std::string& name);
	void TryRelease();
	bool loaded;
	int n;
	Vec3d* points;
	Vec3d* normals;
	Octree* octree;
	Vec3d minPoint, maxPoint;
};

FileLoader::FileLoader() {
	loaded = false;
	n = 0;
	points = NULL;
	normals = NULL;
	octree = NULL;
	minPoint = Vec3d();
	maxPoint = Vec3d();
}

FileLoader::~FileLoader() {
	TryRelease();
}

bool FileLoader::Load(const std::string& name) {
	std::cout << "Loading: " << name << std::endl;
	TryRelease();//初始化
	FILE* fp;
	fopen_s(&fp, ("data/" + name).c_str(), "r");
	if (!fp) { return false; }
	Vec3d temp;
	// 读取数据点数量n
	fscanf_s(fp, "%d", &n);
	points = new Vec3d[n];
	normals = new Vec3d[n];
	// 读取点云并计算点云范围
	fscanf_s(fp, "%lf %lf %lf", &temp.x, &temp.y, &temp.z);
	points[0] = minPoint = maxPoint = temp;
	for (int i = 1; i < n; i++) {
		fscanf_s(fp, "%lf %lf %lf", &temp.x, &temp.y, &temp.z);
		points[i] = temp;
		minPoint.x = (minPoint.x <= temp.x) ? (minPoint.x) : (temp.x);
		minPoint.y = (minPoint.y <= temp.y) ? (minPoint.y) : (temp.y);
		minPoint.z = (minPoint.z <= temp.z) ? (minPoint.z) : (temp.z);
		maxPoint.x = (maxPoint.x >= temp.x) ? (maxPoint.x) : (temp.x);
		maxPoint.y = (maxPoint.y >= temp.y) ? (maxPoint.y) : (temp.y);
		maxPoint.z = (maxPoint.z >= temp.z) ? (maxPoint.z) : (temp.z);
	}
	// 读取法向量
	for (int i = 0; i < n; i++) {
		fscanf_s(fp, "%lf %lf %lf", &temp.x, &temp.y, &temp.z);
		normals[i] = temp;
	}
	fclose(fp);
	// 创建八叉树
	octree = new Octree(minPoint, maxPoint, 8);
	octree->AddPoints(points, n);
	loaded = true;
	printf("minPoint:[%f %f %f]\n", minPoint.x, minPoint.y, minPoint.z);
	printf("maxPoint:[%f %f %f]\n", maxPoint.x, maxPoint.y, maxPoint.z);
	printf("Load finished.\n");
	return true;
}

void FileLoader::TryRelease() {
	loaded = false;
	n = 0;
	if (points != NULL) { delete[] points; }
	points = NULL;
	if (normals != NULL) { delete[] normals; }
	normals = NULL;
	if (octree != NULL) { delete octree; }
	octree = NULL;
}