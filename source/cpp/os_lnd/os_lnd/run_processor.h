#pragma once
#include "lu_run_strategy.h"
#include "odeint_rk4_run_strategy.h"
#include "smallest_eigen_vector_run_strategy.h"
#include "all_evals_run_strategy.h"
#include "eigen_dense_run_strategy.h"
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
		else if (task == "odeint")
		{
			run_strategy = std::make_unique<ODEIntRK4RunStrategy>();
		}
		else if (task == "smallest_eigen_vector")
		{
			run_strategy = std::make_unique<SmallestEigenVectorRunStrategy>();
		}
		else if (task == "all_evals")
		{
			run_strategy = std::make_unique<AllEvalsRunStrategy>();
		}
		else if (task == "eigen_dense")
		{
			run_strategy = std::make_unique<EigenDenseRunStrategy>();
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

	void process_serial(Model& model, std::map<std::string, std::vector<double>>& features_double, std::map<std::string, std::vector<std::complex<double>>>& features_complex)
	{
		
	}
};
