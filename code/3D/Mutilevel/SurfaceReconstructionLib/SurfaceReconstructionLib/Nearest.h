#pragma once
#include <vector>
class Nearest
{
public:
	inline void Deletenearest();
	Nearest();
	double r;
	std::vector<int>* indices;  // 近邻信息
	std::vector<double>* squaredDists;  // 近邻对应距离信息
	int count; //支撑半径内点的个数
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
