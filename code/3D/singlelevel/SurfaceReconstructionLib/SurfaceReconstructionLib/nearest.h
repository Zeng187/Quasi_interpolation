#pragma once
class Nearest
{
public:
	inline void Deletenearest();
	Nearest();
	std::vector<int>* indices;  // ������Ϣ
	std::vector<double>* squaredDists;  // ���ڶ�Ӧ������Ϣ
	int count; //֧�Ű뾶�ڵ�ĸ���
	double r;
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
