#pragma once
#include "Vec3d.h"
#include <queue>
#include <vector>

class CubeNode
{
public:

	std::vector<int>* pointIds;  // 叶结点包含点的序号表的指针，非叶结点为NULL
	Vec3d centroid;
	Vec3d Normal;
	int count; //记录方块内部的点数量
	Vec3d minPoint;
	Vec3d maxPoint;
	Vec3d midPoint;
	CubeNode* childs[8];  // 非叶结点的子结点指针表，不存在的子结点为NULL
	int thislevel;
	bool isleaf;
	bool complete;
	CubeNode* next;//用于构建横向链
	CubeNode();
	CubeNode(Vec3d minPoint, Vec3d maxPoint, int tlevel, int level);
	~CubeNode();
	int ConstructCubeNodelink(CubeNode* &h);
	void Collect(Vec3d* Dpoints, Vec3d* Dnormals, CubeNode* h, int n);//求取层次网格第d层的对应节点和法向量，返回点的数量；
};

CubeNode::CubeNode()
{
	this->next = NULL;
	this->pointIds = NULL;

}

CubeNode::CubeNode(Vec3d minPoint, Vec3d maxPoint,int tlevel,int level)
{
	this->minPoint = minPoint;
	this->maxPoint = maxPoint;
	this->midPoint = (minPoint + maxPoint) / 2;
	this->thislevel = tlevel;
	this->isleaf = tlevel < level ? false : true;
	pointIds = new std::vector<int>;
	for (int i = 0; i < 8; i++) {
		childs[i] = NULL;
	}
	this->count = 0;
	this->complete = false;
	this->next = NULL;
}

CubeNode::~CubeNode()
{
	if (isleaf) {
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

class CubeTree
{
public:
	int level;
	std::vector<Vec3d> points;  // 存储所有点的坐标，排列顺序为点的添加顺序
	Vec3d* normals;
	CubeNode* root;
	CubeTree(Vec3d minPoint, Vec3d maxPoint, Vec3d* normals, int level);
	~CubeTree();
	void AddPoint(Vec3d point);  // 为八叉树添加点
	void AddPoints(Vec3d* points, int n);  // 为八叉树添加点集
	void Caculatecentriodandnormals(CubeNode*&h);
private:
	void AddPointToNode(CubeNode* node, int pointId);  // 为结点添加点
};

CubeTree::CubeTree(Vec3d minPoint, Vec3d maxPoint, Vec3d* normals,int level)
{
	this->level = level;
	root = new CubeNode(minPoint, maxPoint,0,level);
	this->normals = normals;
}

CubeTree::~CubeTree()
{
	delete root;
}


void CubeTree::AddPoint(Vec3d point) {
	points.push_back(point);
	AddPointToNode(root, points.size() - 1);
}

void CubeTree::AddPoints(Vec3d* points, int n) {
	for (int i = 0; i < n; i++) {
		AddPoint(points[i]);
	}
}

void CubeTree::AddPointToNode(CubeNode* node, int pointId)
{
	if (node->isleaf)
	{
		node->pointIds->push_back(pointId);
		node->count++;
	}
	else
	{
		int subNodeId = 0;
		if (points[pointId].x > node->midPoint.x) { subNodeId += 1; }
		if (points[pointId].y > node->midPoint.y) { subNodeId += 2; }
		if (points[pointId].z > node->midPoint.z) { subNodeId += 4; }
		if (node->childs[subNodeId] == NULL) {
			Vec3d subNodeMinPoint = node->minPoint;
			Vec3d subNodeMaxPoint = node->maxPoint;
			if (points[pointId].x > node->midPoint.x) { subNodeMinPoint.x = node->midPoint.x; }
			else { subNodeMaxPoint.x = node->midPoint.x; }
			if (points[pointId].y > node->midPoint.y) { subNodeMinPoint.y = node->midPoint.y; }
			else { subNodeMaxPoint.y = node->midPoint.y; }
			if (points[pointId].z > node->midPoint.z) { subNodeMinPoint.z = node->midPoint.z; }
			else { subNodeMaxPoint.z = node->midPoint.z; }
			node->childs[subNodeId] = new CubeNode(subNodeMinPoint, subNodeMaxPoint, node->thislevel + 1, level);
		}
		AddPointToNode(node->childs[subNodeId], pointId);
	}
	
}

void CubeTree::Caculatecentriodandnormals(CubeNode *&h)
{
	Vec3d c = Vec3d(); Vec3d n = Vec3d(); int i;
	if (!h->complete)
	{
		if (h->isleaf)
		{
			
			for ( i= 0; i < h->count; i++)
			{
				//std::cout << points[(*h->pointIds)[i]].x << std::endl;
				c= c + points[(*h->pointIds)[i]];
				//std::cout << c.x << std::endl;
				n= n + normals[(*h->pointIds)[i]];
			}
			//std::cout << c.x << std::endl;
			//std::cout << c.x/h->count << std::endl;
			h->centroid = c / h->count;
			//std::cout << h->centroid.x << std::endl;
			h->Normal = n / h->count;
			h->complete = true;
		}
		else {
			//当前节点不是叶子节点
			for (i = 0; i < 8; i++)
			{
				if (h->childs[i] != NULL)
				{
					if ((!(h->childs[i])->complete))
						Caculatecentriodandnormals(h->childs[i]);
					h->count = h->count + h->childs[i]->count;
					c = c + h->childs[i]->centroid * h->childs[i]->count;
					n = n + h->childs[i]->Normal * h->childs[i]->count;
				}
			}
			h->centroid = c / h->count;
			h->Normal = n / h->count;
			h->complete = true;
		}

	}
}




int CubeNode::ConstructCubeNodelink(CubeNode* &h)
{
	//由当前CubeNodelink构建下一层CubeNodelink；
	CubeNode* r = h, * p = h; int count = 0;
	while (p->next != NULL)
	{
		p = p->next;
		for (int i = 0; i < 8; i++)
		{
			if (p->childs[i] != NULL)
			{
				count = count + 1;
				r->next = p->childs[i];
				r = r->next;
			}
		}
	}
	r->next = NULL;
	
	return count;
}

void CubeNode::Collect(Vec3d* Dpoints, Vec3d* Dnormals, CubeNode* h,int n)
{

	CubeNode* p = h->next;
	for (int i = 0; i < n; i++)
	{
		Dpoints[i] = p->centroid;
		Dnormals[i] = p->Normal;
		p = p->next;
	}
}



