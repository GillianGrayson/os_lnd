#pragma once
#include "lindbladian_lu_model_run_strategy.h"
#include "lindbladian_odeint_rk4_model_run_strategy.h"


struct ModelRunProcessor
{
	std::unique_ptr<ModelRunStrategy> model_strategy;

	void set_model_strategy(Model& model)
	{
		const std::string system = model.ini.Get("global", "task", "unknown");
		if (system == "lindbladian_lu")
		{
			model_strategy = std::make_unique<LindbladianLUModelRunStrategy>();
		}
		else if (system == "lindbladian_odeint_rk4")
		{
			model_strategy = std::make_unique<LindbladianODEIntRK4ModelRunStrategy>();
		}
		else
		{
			model.throw_error("Unsupported model strategy");
		}
	}

	void process_model(Model& model) const
	{
		model_strategy->run(model);
	}
};
