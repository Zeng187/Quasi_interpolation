// SurfaceReconstructionLib.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include "Octree.h"
#include "Vec3d.h"
#include "Mat3.h"
#include "IRBF.h"
#include "FileLoader.h"
#include "ResultSaver.h"
#include "Batch.h"
#include <cstdio>
#include <string>
using namespace std;

void DeepFirstTravel(OctreeNode* node) {
	if (node->isLeaf) {
		printf("MinPoint:(%f,%f,%f),MaxPoint:(%f,%f,%f),Points:",
			node->minPoint.x, node->minPoint.y, node->minPoint.z,
			node->maxPoint.x, node->maxPoint.y, node->maxPoint.z
		);
		for (int i = 0; i < node->pointIds->size(); i++) {
			printf("%d,", (*(node->pointIds))[i]);
		}
		printf("\n");
	}
	else {
		printf("MinPoint:(%f,%f,%f),MaxPoint:(%f,%f,%f),Not a leaf\n",
			node->minPoint.x, node->minPoint.y, node->minPoint.z,
			node->maxPoint.x, node->maxPoint.y, node->maxPoint.z
		);
		for (int i = 0; i < 8; i++) {
			if (node->childs[i] != NULL) {
				DeepFirstTravel(node->childs[i]);
			}
		}
	}
}

void KNNTest(Octree* octree) {
	const int k = 9;
	int pointIds[k];
	double squaredDistance[k];
	Vec3d targetPoint(0, 0, 0);
	int count = octree->FindKNearestNeighbors(targetPoint, k, pointIds, squaredDistance);
	printf("Find %d results:\n", count);
	for (int i = 0; i < count; i++) {
		printf("%d,%f\n", pointIds[i], squaredDistance[i]);
	}
	int targetid = 8;
	count = octree->FindKNearestNeighbors(targetid, k, pointIds, squaredDistance);
	printf("Find %d results:\n", count);
	for (int i = 0; i < count; i++) {
		printf("%d,%f\n", pointIds[i], squaredDistance[i]);
	}
}

void LinearSolverTest() {
	Mat3 A = Mat3(
		Vec3d(0.001, -1, -2),
		Vec3d(2, 3.712, 1.072),
		Vec3d(3, 4.623, 5.643)
	);
	Vec3d b = Vec3d(1, 2, 3);
	Vec3d x = A.SolveLinerEquation(b);
	printf("%f, %f, %f\n", x.x, x.y, x.z);
}




int main()
{
	printf("Input batch file name: ");
	string name; int storageway;
	std::cin >> name;
	printf("Input the way of storage:(1 is mat, 2 is ply) ");
	std::cin >> storageway;
	Batch* batch = new Batch(name, storageway);
	batch->Execute();
	delete batch;
	system("pause");
}