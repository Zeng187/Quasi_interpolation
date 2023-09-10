#pragma once
#include <vector>
class Nearest
{
public:
	inline void Deletenearest();
	Nearest();
	double r;
	std::vector<int>* indices;  // ������Ϣ
	std::vector<double>* squaredDists;  // ���ڶ�Ӧ������Ϣ
	int count; //֧�Ű뾶�ڵ�ĸ���
};

Nearest::Nearest()
{
	indices = new std::vector<int>();
	squaredDists = new std::vector<double>();
}

inline void Nearest::Deletenearest()
{
	delete indices;
	delete squaredDists;
}
