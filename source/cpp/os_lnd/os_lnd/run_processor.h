#pragma once
#include "lindbladian_lu_run_strategy.h"
#include "lindbladian_odeint_rk4_run_strategy.h"
#include "lindbladian_smaller_eigen_vector_run_strategy.h"


struct RunProcessor
{
	std::unique_ptr<RunStrategy> run_strategy;

	void set_strategy(Model& model)
	{
		const std::string task = model.ini.Get("global", "task", "unknown");
		if (task == "lindbladian_lu")
		{
			run_strategy = std::make_unique<LindbladianLURunStrategy>();
		}
		else if (task == "lindbladian_odeint_rk4")
		{
			run_strategy = std::make_unique<LindbladianODEIntRK4RunStrategy>();
		}
		else if (task == "lindbladian_smaller_eigen_vector")
		{
			run_strategy = std::make_unique<LindbladianSmallerEigenVectorRunStrategy>();
		}
		else
		{
			model.throw_error("Unsupported task");
		}
	}

	void process(Model& model) const
	{
		run_strategy->run(model);
	}
};
