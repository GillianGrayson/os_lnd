#pragma once
#include "init.h"

inline sp_mtx get_sp_eye(const int size)
{
	std::vector<triplet> vec_triplets;
	vec_triplets.reserve(size);
	for (int st_id = 0; st_id < size; st_id++)
	{
		vec_triplets.push_back(triplet(st_id, st_id, std::complex<double>(1.0, 0.0)));
	}

	sp_mtx mtx(size, size);
	mtx.setFromTriplets(vec_triplets.begin(), vec_triplets.end());

	return mtx;
}