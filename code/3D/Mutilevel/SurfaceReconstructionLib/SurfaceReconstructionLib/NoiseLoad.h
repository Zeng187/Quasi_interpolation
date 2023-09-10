#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <cstdio>
#include "Vec3d.h"
#include "Octree.h"
#include "Cube.h"

class NoiseLoader {
public:
	NoiseLoader();
	~NoiseLoader();
	bool Load(const std::string& name);
	void TryRelease();
	bool loaded;
	int n;
	int* noise;
};

NoiseLoader::NoiseLoader() {
	n = 0;
	noise = NULL;
}

NoiseLoader::~NoiseLoader() {
	TryRelease();
}

bool NoiseLoader::Load(const std::string& name) {
	std::cout << "Loading: " << name << std::endl;
	TryRelease();//初始化
	FILE* fp;
	fopen_s(&fp, ("data/" + name).c_str(), "r");
	if (!fp) { return false; }
	int temp;
	// 读取噪声点数量n
	fscanf_s(fp, "%d", &n);
	noise = new int[n];
	// 读取噪声序号
	
	for (int i = 0; i < n; i++) {
		fscanf_s(fp, "%d", &temp);
		noise[i] = temp;
	}
	fclose(fp);
	// 创建八叉树
	printf("Load finished.\n");
	return true;
}

void NoiseLoader::TryRelease() {
	n = 0;
	if (noise!= NULL) { delete[] noise; }
	noise = NULL;
}