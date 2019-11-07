#pragma once
#include "model_run_strategy.h"
#include <Eigen/SparseLU>
#include "save.h"

struct LUModelRunStrategy : ModelRunStrategy
{
	void run_asymptotic_rho(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);

		Eigen::SparseLU<sp_mtx> solver;
		solver.compute(model.lindbladian);
		if (solver.info() != Eigen::Success)
		{
			throw std::runtime_error("Decomposition failed");
		}
		else
		{
			std::cout << "Decomposition complete!" << std::endl;
		}

		Eigen::VectorXcd right_part = Eigen::VectorXcd::Zero(model.sys_size * model.sys_size);
		Eigen::VectorXcd rho_vec = solver.solve(right_part);
		if (solver.info() != Eigen::Success)
		{
			throw std::runtime_error("Solving failed");
		}
		else
		{
			std::cout << "Solving complete!" << std::endl;
		}

		model.rho = Eigen::Map<Eigen::MatrixXcd>(rho_vec.data(), model.sys_size, model.sys_size);

		auto fn = "rho_mtx" + model.suffix;
		save_dense(model.rho, fn, save_precision);
	}
};