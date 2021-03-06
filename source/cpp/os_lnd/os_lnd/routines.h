#pragma once
#include "init.h"
#include <unsupported/Eigen/KroneckerProduct>

inline sp_mtx get_sigma_x()
{
	std::vector<triplet> vec_triplets;
	vec_triplets.push_back(triplet(0, 1, { 1.0, 0.0 }));
	vec_triplets.push_back(triplet(1, 0, { 1.0, 0.0 }));
	sp_mtx sigma(2, 2);
	sigma.setFromTriplets(vec_triplets.begin(), vec_triplets.end());

	return sigma;
}

inline sp_mtx get_sigma_y()
{
	std::vector<triplet> vec_triplets;
	vec_triplets.push_back(triplet(0, 1, { 0.0, -1.0 }));
	vec_triplets.push_back(triplet(1, 0, { 0.0, 1.0 }));
	sp_mtx sigma(2, 2);
	sigma.setFromTriplets(vec_triplets.begin(), vec_triplets.end());

	return sigma;
}

inline sp_mtx get_sigma_z()
{
	std::vector<triplet> vec_triplets;
	vec_triplets.push_back(triplet(0, 0, { 1.0, 0.0 }));
	vec_triplets.push_back(triplet(1, 1, { -1.0, 0.0 }));
	sp_mtx sigma(2, 2);
	sigma.setFromTriplets(vec_triplets.begin(), vec_triplets.end());

	return sigma;
}

inline sp_mtx get_sigma_m()
{
	std::vector<triplet> vec_triplets;
	vec_triplets.push_back(triplet(1, 0, { 1.0, 0.0 }));
	sp_mtx sigma(2, 2);
	sigma.setFromTriplets(vec_triplets.begin(), vec_triplets.end());

	return sigma;
}

inline sp_mtx get_sigma_p()
{
	std::vector<triplet> vec_triplets;
	vec_triplets.push_back(triplet(0, 1, { 1.0, 0.0 }));
	sp_mtx sigma(2, 2);
	sigma.setFromTriplets(vec_triplets.begin(), vec_triplets.end());

	return sigma;
}

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

inline std::vector<sp_mtx> get_kronecker_mtxs(const int num_particles, std::string mtx_type)
{
	std::vector<sp_mtx> kronecker_mtxs;

	sp_mtx sigma_0 = get_sp_eye(2);
	sp_mtx sigma_x = get_sigma_x();
	sp_mtx sigma_y = get_sigma_y();
	sp_mtx sigma_z = get_sigma_z();
	sp_mtx sigma_p = get_sigma_p();
	sp_mtx sigma_m = get_sigma_m();

	sp_mtx target_mtx;
	if (mtx_type == "sigma_0")
	{
		target_mtx = sigma_0;
	}
	else if (mtx_type == "sigma_x")
	{
		target_mtx = sigma_x;
	}
	else if (mtx_type == "sigma_y")
	{
		target_mtx = sigma_y;
	}
	else if (mtx_type == "sigma_z")
	{
		target_mtx = sigma_z;
	}
	else if (mtx_type == "sigma_p")
	{
		target_mtx = sigma_p;
	}
	else if (mtx_type == "sigma_m")
	{
		target_mtx = sigma_m;
	}
	else
	{
		throw std::runtime_error("Unsupported kronecker matrix type");
	}


	for (auto spin_id = 0; spin_id < num_particles; spin_id++)
	{
		sp_mtx curr_mtx;

		if (spin_id == 0)
		{
			curr_mtx = target_mtx;
		}
		else
		{
			curr_mtx = sigma_0;
		}

		for (auto inner_id = 1; inner_id < num_particles; inner_id++)
		{
			if (inner_id == spin_id)
			{
				curr_mtx = Eigen::kroneckerProduct(curr_mtx, target_mtx).eval();
			}
			else
			{
				curr_mtx = Eigen::kroneckerProduct(curr_mtx, sigma_0).eval();
			}
		}
		kronecker_mtxs.push_back(curr_mtx);
	}

	if (kronecker_mtxs.size() != num_particles)
	{
		throw std::runtime_error("Error during kronecker matrices initialization");
	}
	
	return kronecker_mtxs;
}

inline ds_mtx get_reshuffle_ds_mtx_1(ds_mtx& mtx, int full_size, int part_size)
{
	if (full_size != part_size * part_size)
	{
		throw std::runtime_error("reshuffle: full_size must be equal to square of part_size");
	}
	if (mtx.rows() != full_size || mtx.cols() != full_size)
	{
		throw std::runtime_error("reshuffle: matrix is not square");
	}

	ds_mtx reshuffle = ds_mtx::Zero(full_size, full_size);

	for (auto i = 0; i < part_size; i++)
	{
		for (auto j = 0; j < part_size; j++)
		{
			for (auto k = 0; k < part_size; k++)
			{
				for (auto m = 0; m < part_size; m++)
				{
					auto origin_row = i * part_size + k;
					auto origin_col = j * part_size + m;

					auto target_row = i * part_size + j;
					auto target_col = k * part_size + m;

					reshuffle(target_row, target_col) = mtx(origin_row, origin_col);
				}
			}
		}
	}

	return reshuffle;
}

inline ds_mtx get_reshuffle_ds_mtx_0(ds_mtx& mtx, int full_size, int part_size)
{
	if (full_size != part_size * part_size)
	{
		throw std::runtime_error("reshuffle: full_size must be equal to square of part_size");
	}
	if (mtx.rows() != full_size || mtx.cols() != full_size)
	{
		throw std::runtime_error("reshuffle: matrix is not square");
	}

	ds_mtx reshuffle = ds_mtx::Zero(full_size, full_size);

	for (auto i = 0; i < part_size; i++)
	{
		for (auto j = 0; j < part_size; j++)
		{
			for (auto k = 0; k < part_size; k++)
			{
				for (auto m = 0; m < part_size; m++)
				{
					auto origin_row = i * part_size + k;
					auto origin_col = j * part_size + m;

					auto target_row = m * part_size + k;
					auto target_col = j * part_size + i;

					reshuffle(target_row, target_col) = mtx(origin_row, origin_col);
				}
			}
		}
	}

	return reshuffle;
}

inline ds_mtx get_addition_ds_mtx(ds_mtx& mtx, int full_size, int part_size)
{
	if (full_size != part_size * part_size)
	{
		throw std::runtime_error("addition: full_size must be equal to square of part_size");
	}
	if (mtx.rows() != full_size || mtx.cols() != full_size)
	{
		throw std::runtime_error("addition: matrix is not square");
	}

	ds_mtx addition = ds_mtx::Zero(part_size, part_size);

	for (auto s1 = 0; s1 < part_size; s1++)
	{
		for (auto s2 = 0; s2 < part_size; s2++)
		{
			for (auto s3 = 0; s3 < part_size; s3++)
			{
				auto w1 = s3 + part_size * s1;
				auto w2 = s3 + part_size * s2;

				addition(s1, s2) += mtx(w1, w2);
			}
		}
	}

	return addition;
}