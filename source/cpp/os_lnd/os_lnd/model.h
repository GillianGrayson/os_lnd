#pragma once
#include <INIReader.h>
#include <Eigen/SparseCore>
#include <vector>
#include "init.h"

struct Model
{
	INIReader ini;
	std::string suffix;
	int sys_size;
	sp_mtx hamiltonian;
	std::vector<sp_mtx> dissipators;
	sp_mtx lindbladian;

	Model(INIReader& ini) : ini(ini)
	{
	}
};
