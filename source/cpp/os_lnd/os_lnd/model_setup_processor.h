#pragma once
#include "mbl_model_setup_strategy.h"


struct ModelSetupProcessor
{
	std::unique_ptr<ModelSetupStrategy> model_strategy;

	void set_model_strategy(Model& model)
	{
		const std::string system = model.ini.Get("global", "system", "unknown");
		if (system == "mbl")
		{
			model_strategy = std::make_unique<MBLModelSetupStrategy>();
		}
		else
		{
			throw std::runtime_error("Unsupported model strategy");
		}
	}

	void process_model(Model& model) const
	{
		model_strategy->setup_suffix(model);
		model_strategy->setup_sys_size(model);
		model_strategy->setup_hamiltonian(model);
		model_strategy->setup_dissipators(model);
		model_strategy->setup_lindbladian(model);
	}
};