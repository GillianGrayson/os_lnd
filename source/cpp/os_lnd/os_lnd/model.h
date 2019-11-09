#pragma once
#include <INIReader.h>
#include <vector>
#include "init.h"
#include <chrono>

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
	std::chrono::high_resolution_clock::time_point run_time;

	Model(INIReader& ini) : ini(ini)
	{
		run_time = std::chrono::high_resolution_clock::now();
	}
};
