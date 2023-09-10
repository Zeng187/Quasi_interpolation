#pragma once
#include "Vec3d.h"
#include <vector>

struct OctreeNode {
public:
	Vec3d minPoint, maxPoint, midPoint;  // 体块范围（左闭右开）
	bool isLeaf;  // 是否为叶结点
	std::vector<int>* pointIds;  // 叶结点包含点的序号表的指针，非叶结点为NULL
	OctreeNode* childs[8];  // 非叶结点的子结点指针表，不存在的子结点为NULL
	int count;
	OctreeNode(Vec3d minPoint, Vec3d maxPoint);
	~OctreeNode();
	double SquaredDistanceToPoint(Vec3d point);  // 求给定点到体块的最近距离
};

class Octree {
public:
	Octree(Vec3d minPoint, Vec3d maxPoint, int maxPointCountPerNode);  // 根据给定范围（左闭右开）创建八叉树
	~Octree();
	std::vector<Vec3d> points;  // 存储所有点的坐标，排列顺序为点的添加顺序
	void AddPoint(Vec3d point);  // 为八叉树添加点
	void AddPoints(Vec3d* points, int n);  // 为八叉树添加点集
	int FindKNearestNeighbors(int pointId, int k, int* pointIds, double* squaredDistance); // 寻找八叉树内点的k近邻，排列顺序由小到大
	int FindKNearestNeighbors(Vec3d targetPoint, int k, int* pointIds, double* squaredDistance); // 寻找任意点的k近邻，排列顺序由小到大
	int FindNearestNeighborsBySquaredRadius(Vec3d targetPoint, double squaredRadius, std::vector<int>* pointIds, std::vector<double>* squaredDistance, double& averagedistance); // 寻找任意点的k近邻，排列顺序由小到大
	int FindNearestNeighborsBySquaredRadius(Vec3d targetPoint, double squaredRadius, std::vector<int>* pointIds, std::vector<double>* squaredDistance); // 寻找任意点的k近邻，排列顺序由小到大
	OctreeNode* root;  // 八叉树根节点
	double Averageleaflength(OctreeNode* r, int &count);
private:
	int maxPointCountPerNode; // 体块包含点最大个数
	void AddPointToNode(OctreeNode* node, int pointId);  // 为结点添加点
	void SplitNode(OctreeNode* node); // 分裂结点
	
};

struct OctreeSearchNode {
	bool isPoint;  // 该搜索结点是否为一个点
	int pointId;  // 该搜索结点对应的坐标点
	OctreeNode* node;  // 该搜索结点对应的八叉树体块
	double squaredDistanceToTarget;  // 该搜索结点与目标点间的平方距离
	inline bool operator < (const OctreeSearchNode& node) const {
		return squaredDistanceToTarget > node.squaredDistanceToTarget;
	}
	OctreeSearchNode(OctreeNode* node, double squaredDistance);
	OctreeSearchNode(int pointId, double squaredDistance);
};