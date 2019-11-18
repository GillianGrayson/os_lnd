#pragma once
#include "run_strategy.h"
#include <Eigen/SparseLU>
#include "save.h"

struct LindbladianLURunStrategy : RunStrategy
{
	void run(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);

		Eigen::SparseLU<sp_mtx> solver;
		solver.compute(model.lindbladian);
		if (solver.info() != Eigen::Success)
		{
			model.throw_error("Decomposition failed");
		}
		else
		{
			model.log_message("Decomposition complete!");
		}

		Eigen::VectorXcd right_part = Eigen::VectorXcd::Zero(model.sys_size * model.sys_size);
		Eigen::VectorXcd rho_vec = solver.solve(right_part);
		if (solver.info() != Eigen::Success)
		{
			model.throw_error("Solving failed");
		}
		else
		{
			model.log_message("Solving complete!");
		}

		model.rho = Eigen::Map<Eigen::MatrixXcd>(rho_vec.data(), model.sys_size, model.sys_size);
		auto fn = "rho_mtx" + model.suffix;
		save_dense_mtx(model.rho, fn, save_precision);
	}
};