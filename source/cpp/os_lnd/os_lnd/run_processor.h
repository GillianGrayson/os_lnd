#pragma once
#include "lindbladian_lu_run_strategy.h"
#include "lindbladian_odeint_rk4_run_strategy.h"


struct ModelRunProcessor
{
	std::unique_ptr<RunStrategy> run_strategy;

	void set_model_strategy(Model& model)
	{
		const std::string system = model.ini.Get("global", "task", "unknown");
		if (system == "lindbladian_lu")
		{
			run_strategy = std::make_unique<LindbladianLURunStrategy>();
		}
		else if (system == "lindbladian_odeint_rk4")
		{
			run_strategy = std::make_unique<LindbladianODEIntRK4RunStrategy>();
		}
		else
		{
			model.throw_error("Unsupported run_strategy");
		}
	}

	void process_model(Model& model) const
	{
		run_strategy->run(model);
	}
};
