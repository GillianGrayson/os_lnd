#pragma once
#include "mbl_strategy.h"


struct ModelProcessor
{
	std::unique_ptr<ModelStrategy> model_strategy;

	void set_model_strategy(Model& model)
	{
		const std::string system = model.ini.Get("global", "system", "unknown");
		if (system == "mbl")
		{
			model_strategy = std::make_unique<MBLModelStrategy>();
		}
		else
		{
			throw std::runtime_error("Unsupported model strategy");
		}
	}

	void process_model(Model& model) const
	{
		model_strategy->set_suffix(model);
		model_strategy->set_sys_size(model);
		model_strategy->set_hamiltonian(model);
	}
};