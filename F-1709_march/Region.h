#pragma once
#include <vector>

struct Region
{
public:
	int N;
	int material;
	std::vector<int> def_nodes;
	std::vector<int> inside_nodes_B;
	std::vector<int> inside_nodes_E;
};