#include "pch.h"
#include "Octree.h"
#include "Vec3d.h"
#include <vector>
#include <queue>

OctreeNode::OctreeNode(Vec3d minPoint, Vec3d maxPoint) {
	this->minPoint = minPoint;
	this->maxPoint = maxPoint;
	midPoint = (minPoint + maxPoint) / 2;
	isLeaf = true;
	pointIds = new std::vector<int>;
	for (int i = 0; i < 8; i++) {
		childs[i] = NULL;
	}
}
OctreeNode::~OctreeNode() {
	if (isLeaf) {
		delete pointIds;
	}
	else {
		for (int i = 0; i < 8; i++) {
			if (childs[i] != NULL) {
				delete childs[i];
			}
		}
	}
}
double OctreeNode::SquaredDistanceToPoint(Vec3d point) {
	Vec3d differVec;  // 点到体块上最近点的距离向量
	if (point.x < minPoint.x) { differVec.x = minPoint.x - point.x; }
	else if (point.x > maxPoint.x){ differVec.x = point.x - maxPoint.x; }
	else { differVec.x = 0; }
	if (point.y < minPoint.y) { differVec.y = minPoint.y - point.y; }
	else if (point.y > maxPoint.y) { differVec.y = point.y - maxPoint.y; }
	else { differVec.y = 0; }
	if (point.z < minPoint.z) { differVec.z = minPoint.z - point.z; }
	else if (point.z > maxPoint.z) { differVec.z = point.z - maxPoint.z; }
	else { differVec.z = 0; }
	return differVec.SquaredLength();
}

Octree::Octree(Vec3d minPoint, Vec3d maxPoint, int maxPointCountPerNode) {
	root = new OctreeNode(minPoint, maxPoint);
	this->maxPointCountPerNode = maxPointCountPerNode;
}
Octree::~Octree() {
	delete root;
}
void Octree::AddPoint(Vec3d point) {
	points.push_back(point);
	AddPointToNode(root, points.size() - 1);
}
void Octree::AddPoints(Vec3d* points, int n) {
	for (int i = 0; i < n; i++) {
		AddPoint(points[i]);
	}
}
void Octree::AddPointToNode(OctreeNode* node, int pointId) {
	if (node->isLeaf) {
		node->pointIds->push_back(pointId);
		if (node->pointIds->size() > maxPointCountPerNode) {  // 体块包含点数量高于最大数量时进行分裂
			SplitNode(node);
		}
	}
	else {
		// 先计算坐标点对应的子体块id，id由0开始，顺序为xyz轴由小到大
		int subNodeId = 0;
		if (points[pointId].x > node->midPoint.x) { subNodeId += 1; }
		if (points[pointId].y > node->midPoint.y) { subNodeId += 2; }
		if (points[pointId].z > node->midPoint.z) { subNodeId += 4; }
		// 对应子体块不存在时先创建
		if (node->childs[subNodeId] == NULL) {
			Vec3d subNodeMinPoint = node->minPoint;
			Vec3d subNodeMaxPoint = node->maxPoint;
			if (points[pointId].x > node->midPoint.x) { subNodeMinPoint.x = node->midPoint.x; }
			else { subNodeMaxPoint.x = node->midPoint.x; }
			if (points[pointId].y > node->midPoint.y) { subNodeMinPoint.y = node->midPoint.y; }
			else { subNodeMaxPoint.y = node->midPoint.y; }
			if (points[pointId].z > node->midPoint.z) { subNodeMinPoint.z = node->midPoint.z; }
			else { subNodeMaxPoint.z = node->midPoint.z; }
			node->childs[subNodeId] = new OctreeNode(subNodeMinPoint, subNodeMaxPoint);
		}
		AddPointToNode(node->childs[subNodeId], pointId);  // 向子体块递归添加坐标点
	}
}
void Octree::SplitNode(OctreeNode* node) {
	node->isLeaf = false;  // 设置结点为非叶结点
	for (int i = 0; i < node->pointIds->size(); i++) {
		AddPointToNode(node, (*(node->pointIds))[i]);
	}
	delete node->pointIds;  // 释放体块坐标点表的空间
}
int Octree::FindKNearestNeighbors(int pointId, int k, int* pointIds, double* squaredDistance) {
	int count = 0;
	Vec3d targetPoint = points[pointId];
	std::priority_queue<OctreeSearchNode> heap;  // 利用优先队列进行k近邻搜索
	heap.push(OctreeSearchNode(root, 0));  // 入队八叉树根节点
	while (!heap.empty() && count < k) {
		OctreeSearchNode node = heap.top();
		heap.pop();
		if (node.isPoint) {  // 当前最近搜索结果为一点，加入结果内
			if (node.pointId == pointId) {  // 目标点自身不加入结果
				continue;
			}
			pointIds[count] = node.pointId;
			squaredDistance[count] = node.squaredDistanceToTarget;
			count++;
		}
		else {  // 当前最近搜索结果为一体块，拆解并加入队列
			if (node.node->isLeaf) {  // 体块为叶结点时，将所含点加入队列
				for (int i = 0; i < node.node->pointIds->size(); i++) {
					int id = (*(node.node->pointIds))[i];
					double squaredDistance = (points[id] - targetPoint).SquaredLength();
					heap.push(OctreeSearchNode(id, squaredDistance));
				}
			}
			else {  // 体块为非叶结点时，将所含点加入队列
				for (int i = 0; i < 8; i++) {
					if (node.node->childs[i] == NULL) { continue; }
					OctreeNode* subNode = node.node->childs[i];
					double squaredDistance = subNode->SquaredDistanceToPoint(targetPoint);
					heap.push(OctreeSearchNode(subNode, squaredDistance));
				}
			}
		}
	}
	return count;
}
int Octree::FindKNearestNeighbors(Vec3d targetPoint, int k, int* pointIds, double* squaredDistance) {
	int count = 0;
	std::priority_queue<OctreeSearchNode> heap;  // 利用优先队列进行k近邻搜索
	heap.push(OctreeSearchNode(root, 0));  // 入队八叉树根节点
	while (!heap.empty() && count < k) {
		OctreeSearchNode node = heap.top();
		heap.pop();
		if (node.isPoint) {  // 当前最近搜索结果为一点，加入结果内
			pointIds[count] = node.pointId;
			squaredDistance[count] = node.squaredDistanceToTarget;
			count++;
		}
		else {  // 当前最近搜索结果为一体块，拆解并加入队列
			if (node.node->isLeaf) {  // 体块为叶结点时，将所含点加入队列
				for (int i = 0; i < node.node->pointIds->size(); i++) {
					int id = (*(node.node->pointIds))[i];
					double squaredDistance = (points[id] - targetPoint).SquaredLength();
					heap.push(OctreeSearchNode(id, squaredDistance));
				}
			}
			else {  // 体块为非叶结点时，将所含点加入队列
				for (int i = 0; i < 8; i++) {
					if (node.node->childs[i] == NULL) { continue; }
					OctreeNode* subNode = node.node->childs[i];
					double squaredDistance = subNode->SquaredDistanceToPoint(targetPoint);
					heap.push(OctreeSearchNode(subNode, squaredDistance));
				}
			}
		}
	}
	return count;
}

int Octree::FindNearestNeighborsBySquaredRadius(Vec3d targetPoint, double squaredRadius, std::vector<int>* pointIds, std::vector<double>* squaredDistance) {
	int count = 0;
	pointIds->clear();
	squaredDistance->clear();
	std::priority_queue<OctreeSearchNode> heap;  // 利用优先队列进行k近邻搜索
	heap.push(OctreeSearchNode(root, 0));  // 入队八叉树根节点
	while (!heap.empty()) {
		OctreeSearchNode node = heap.top();
		heap.pop();
		if (node.isPoint) {  // 当前最近搜索结果为一点，加入结果内
			pointIds->push_back(node.pointId);
			squaredDistance->push_back(node.squaredDistanceToTarget);
			count++;
		}
		else {  // 当前最近搜索结果为一体块，拆解并将所求半径内节点加入队列
			if (node.node->isLeaf) {  // 体块为叶结点时，将所含点加入队列
				for (int i = 0; i < node.node->pointIds->size(); i++) {
					int id = (*(node.node->pointIds))[i];
					double squaredDistance = (points[id] - targetPoint).SquaredLength();
					if (squaredDistance > squaredRadius) { continue; }
					heap.push(OctreeSearchNode(id, squaredDistance));
				}
			}
			else {  // 体块为非叶结点时，将所含点加入队列
				for (int i = 0; i < 8; i++) {
					if (node.node->childs[i] == NULL) { continue; }
					OctreeNode* subNode = node.node->childs[i];
					double squaredDistance = subNode->SquaredDistanceToPoint(targetPoint);
					if (squaredDistance > squaredRadius) { continue; }
					heap.push(OctreeSearchNode(subNode, squaredDistance));
				}
			}
		}
	}
	return count;
}

int Octree::FindNearestNeighborsBySquaredRadius(Vec3d targetPoint, double squaredRadius, std::vector<int>* pointIds, std::vector<double>* squaredDistance, double& averagedistance) {
	int count = 0;
	averagedistance = 0;
	pointIds->clear();
	squaredDistance->clear();
	std::priority_queue<OctreeSearchNode> heap;  // 利用优先队列进行k近邻搜索
	heap.push(OctreeSearchNode(root, 0));  // 入队八叉树根节点
	while (!heap.empty()) {
		OctreeSearchNode node = heap.top();
		heap.pop();
		if (node.isPoint) {  // 当前最近搜索结果为一点，加入结果内
			pointIds->push_back(node.pointId);
			squaredDistance->push_back(node.squaredDistanceToTarget);
			averagedistance += node.squaredDistanceToTarget;
			count++;
		}
		else {  // 当前最近搜索结果为一体块，拆解并将所求半径内节点加入队列
			if (node.node->isLeaf) {  // 体块为叶结点时，将所含点加入队列
				for (int i = 0; i < node.node->pointIds->size(); i++) {
					int id = (*(node.node->pointIds))[i];
					double squaredDistance = (points[id] - targetPoint).SquaredLength();
					if (squaredDistance > squaredRadius) { continue; }
					heap.push(OctreeSearchNode(id, squaredDistance));
				}
			}
			else {  // 体块为非叶结点时，将所含点加入队列
				for (int i = 0; i < 8; i++) {
					if (node.node->childs[i] == NULL) { continue; }
					OctreeNode* subNode = node.node->childs[i];
					double squaredDistance = subNode->SquaredDistanceToPoint(targetPoint);
					if (squaredDistance > squaredRadius) { continue; }
					heap.push(OctreeSearchNode(subNode, squaredDistance));
				}
			}
		}
	}
	averagedistance /= count;
	return count;
}

double Octree::Averageleaflength(OctreeNode* r,int& count)
{
	OctreeNode* p = r; double t_result = 0; 
	if (p->isLeaf != true)
	{
		for (int i = 0; i < 8; i++)
		{
			if (p->childs[i] != NULL)
			{
				t_result += Averageleaflength(p->childs[i],count);
			}
				
		}
	}
	else
	{
		count = count + 1;
		return (p->maxPoint - p->minPoint).Length();
	}
	return t_result;
}



OctreeSearchNode::OctreeSearchNode(OctreeNode* node, double squaredDistance) {
	isPoint = false;
	this->node = node;
	squaredDistanceToTarget = squaredDistance;
}
OctreeSearchNode::OctreeSearchNode(int pointId, double squaredDistance) {
	isPoint = true;
	this->pointId = pointId;
	squaredDistanceToTarget = squaredDistance;
}

