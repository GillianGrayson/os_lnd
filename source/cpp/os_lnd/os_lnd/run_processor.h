#pragma once
#include "lu_run_strategy.h"
#include "odeint_rk4_run_strategy.h"
#include "smallest_eigen_vector_run_strategy.h"
#include "mock_run_strategy.h"

struct RunProcessor
{
	std::unique_ptr<RunStrategy> run_strategy;

	void set_strategy(Model& model)
	{
		const std::string task = model.ini.Get("global", "task", "unknown");
		if (task == "lu")
		{
			run_strategy = std::make_unique<LURunStrategy>();
		}
		else if (task == "odeint_rk4")
		{
			run_strategy = std::make_unique<ODEIntRK4RunStrategy>();
		}
		else if (task == "smallest_eigen_vector")
		{
			run_strategy = std::make_unique<SmallestEigenVectorRunStrategy>();
		}
		else if (task == "mock")
		{
			run_strategy = std::make_unique<MockRunStrategy>();
		}
		else
		{
			model.throw_error("Unsupported task");
		}
	}

	void process(Model& model) const
	{
		run_strategy->run(model);
		model.log_memory_usage();
	}
};
