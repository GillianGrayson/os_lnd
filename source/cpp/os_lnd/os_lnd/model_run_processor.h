#pragma once
#include "slss_model_run_strategy.h"


struct ModelRunProcessor
{
	std::unique_ptr<ModelRunStrategy> model_strategy;

	void set_model_strategy(Model& model)
	{
		const std::string system = model.ini.Get("global", "task", "unknown");
		if (system == "lu")
		{
			model_strategy = std::make_unique<LUModelRunStrategy>();
		}
		else
		{
			throw std::runtime_error("Unsupported model strategy");
		}
	}

	void process_model(Model& model) const
	{
		model_strategy->run_asymptotic_rho(model);
	}
};