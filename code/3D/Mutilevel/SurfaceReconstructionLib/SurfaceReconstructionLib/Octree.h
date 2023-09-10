#pragma once
#include "Vec3d.h"
#include <vector>

struct OctreeNode {
public:
	Vec3d minPoint, maxPoint, midPoint;  // ��鷶Χ������ҿ���
	bool isLeaf;  // �Ƿ�ΪҶ���
	std::vector<int>* pointIds;  // Ҷ�����������ű��ָ�룬��Ҷ���ΪNULL
	OctreeNode* childs[8];  // ��Ҷ�����ӽ��ָ��������ڵ��ӽ��ΪNULL
	int count;
	OctreeNode(Vec3d minPoint, Vec3d maxPoint);
	~OctreeNode();
	double SquaredDistanceToPoint(Vec3d point);  // ������㵽�����������
};

class Octree {
public:
	Octree(Vec3d minPoint, Vec3d maxPoint, int maxPointCountPerNode);  // ���ݸ�����Χ������ҿ��������˲���
	~Octree();
	std::vector<Vec3d> points;  // �洢���е�����꣬����˳��Ϊ������˳��
	void AddPoint(Vec3d point);  // Ϊ�˲�����ӵ�
	void AddPoints(Vec3d* points, int n);  // Ϊ�˲�����ӵ㼯
	int FindKNearestNeighbors(int pointId, int k, int* pointIds, double* squaredDistance); // Ѱ�Ұ˲����ڵ��k���ڣ�����˳����С����
	int FindKNearestNeighbors(Vec3d targetPoint, int k, int* pointIds, double* squaredDistance); // Ѱ��������k���ڣ�����˳����С����
	int FindNearestNeighborsBySquaredRadius(Vec3d targetPoint, double squaredRadius, std::vector<int>* pointIds, std::vector<double>* squaredDistance, double& averagedistance); // Ѱ��������k���ڣ�����˳����С����
	int FindNearestNeighborsBySquaredRadius(Vec3d targetPoint, double squaredRadius, std::vector<int>* pointIds, std::vector<double>* squaredDistance); // Ѱ��������k���ڣ�����˳����С����
	OctreeNode* root;  // �˲������ڵ�
	double Averageleaflength(OctreeNode* r, int &count);
private:
	int maxPointCountPerNode; // ��������������
	void AddPointToNode(OctreeNode* node, int pointId);  // Ϊ�����ӵ�
	void SplitNode(OctreeNode* node); // ���ѽ��
	
};

struct OctreeSearchNode {
	bool isPoint;  // ����������Ƿ�Ϊһ����
	int pointId;  // ����������Ӧ�������
	OctreeNode* node;  // ����������Ӧ�İ˲������
	double squaredDistanceToTarget;  // �����������Ŀ�����ƽ������
	inline bool operator < (const OctreeSearchNode& node) const {
		return squaredDistanceToTarget > node.squaredDistanceToTarget;
	}
	OctreeSearchNode(OctreeNode* node, double squaredDistance);
	OctreeSearchNode(int pointId, double squaredDistance);
};