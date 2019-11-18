#pragma once
#include <vector>
#include <fstream>
#include "third_party/eigen-git-mirror/Eigen/Core"

template <class T>
void load_vector(std::vector<T>& v, const std::string& fn)
{
	std::ifstream is(fn);
	T tmp;
	while (is >> tmp)
	{
		v.push_back(tmp);
	}
	is.close();
}
