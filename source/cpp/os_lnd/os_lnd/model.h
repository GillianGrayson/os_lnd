#pragma once
#include <INIReader.h>
#include <vector>
#include "init.h"

struct Model
{
	INIReader ini;
	std::string suffix;
	int sys_size;
	sp_mtx hamiltonian;
	sp_mtx hamiltonian_drv;
	std::vector<sp_mtx> dissipators;
	sp_mtx lindbladian;
	sp_mtx lindbladian_drv;
	Eigen::MatrixXcd rho;

	Model(INIReader& ini) : ini(ini)
	{
	}
};
