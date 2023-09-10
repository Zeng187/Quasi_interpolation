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

/*void IRBFTest() {
	Vec3d points[6] = {
		Vec3d(0, 0, 1),
		Vec3d(0, 0, -1),
		Vec3d(-1, 0, 0),
		Vec3d(0, -1, 0),
		Vec3d(0, 1, 0),
		Vec3d(1, 0, 0),
	};
	Vec3d normals[6] = {
		Vec3d(0, 0, 1),
		Vec3d(0, 0, -1),
		Vec3d(-1, 0, 0),
		Vec3d(0, -1, 0),
		Vec3d(0, 1, 0),
		Vec3d(1, 0, 0),
	};
	Vec3d minPoint(-1, -1, -1);
	Vec3d maxPoint(1, 1, 1);
	Octree* octree = new Octree(minPoint, maxPoint, 16);
	octree->AddPoints(points, 6);
	IRBF* iRBF = new IRBF(5, -1.5, 6, points, normals, octree);
	int res = 6;
	iRBF->SurfaceReconstruction(res, minPoint, maxPoint);
	int index = 0;
	for (int iz = 0; iz <= res; iz++) {
		printf("%d layer:\n", iz);
		for (int ix = 0; ix <= res; ix++) {
			for (int iy = 0; iy <= res; iy++) {
				printf("%f ", iRBF->resultValues[index]);
				index++;
			}
			printf("\n");
		}
	}
	delete iRBF;
	delete octree;
}

void IRBFTest2() {
	string path = "data/bunny453.txt";
	FileLoader* loader = new FileLoader();
	ResultSaver* saver = new ResultSaver();
	loader->Load(path);
	IRBF* iRBF = new IRBF(16, -1.5, loader->n, loader->points, loader->normals, loader->octree);
	int res = 200;
	iRBF->SurfaceReconstruction(res, loader->octree->root->minPoint, loader->octree->root->maxPoint);
	saver->Save(res, iRBF->resultValues, "result");
	delete iRBF;
	delete loader;
	delete saver;
}*/

int main()
{
	printf("Input batch file name: ");
	string name; int storageway;
	std::cin >> name;
	printf("Input the way of storage:(1 is mat, 2 is ply) ");
	std::cin >> storageway;
	Batch* batch = new Batch(name,storageway);
	batch->Execute();
	delete batch;
	system("pause");
}