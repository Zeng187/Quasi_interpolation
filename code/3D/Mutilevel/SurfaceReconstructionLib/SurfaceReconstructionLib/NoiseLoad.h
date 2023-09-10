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
	TryRelease();//��ʼ��
	FILE* fp;
	fopen_s(&fp, ("data/" + name).c_str(), "r");
	if (!fp) { return false; }
	int temp;
	// ��ȡ����������n
	fscanf_s(fp, "%d", &n);
	noise = new int[n];
	// ��ȡ�������
	
	for (int i = 0; i < n; i++) {
		fscanf_s(fp, "%d", &temp);
		noise[i] = temp;
	}
	fclose(fp);
	// �����˲���
	printf("Load finished.\n");
	return true;
}

void NoiseLoader::TryRelease() {
	n = 0;
	if (noise!= NULL) { delete[] noise; }
	noise = NULL;
}