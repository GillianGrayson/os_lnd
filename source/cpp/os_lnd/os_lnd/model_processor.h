#pragma once
#include "mbl_model_strategy.h"
#include "dimer_model_strategy.h"
#include "super_decoh_model_strategy.h"

struct ModelProcessor
{
	std::unique_ptr<ModelStrategy> model_strategy;

	void set_strategy(Model& model)
	{
		const std::string system = model.ini.Get("global", "system", "unknown");
		if (system == "mbl")
		{
			model_strategy = std::make_unique<MBLModelStrategy>();
		}
		else if (system == "dimer")
		{
			model_strategy = std::make_unique<DimerModelStrategy>();
		}
		else if (system == "super_decoh")
		{
			model_strategy = std::make_unique<SuperDecohModelStrategy>();
		}
		else
		{
			model.throw_error("Unsupported model_strategy");
		}
	}

	void init_model(Model& model) const
	{
		model_strategy->setup_suffix(model);
		model_strategy->setup_sys_size(model);
		model_strategy->setup_aux_data(model);
	}

	void release_observables(Model& model) const
	{
		model_strategy->release_observables(model);
	}
	
	void create_model(Model& model) const
	{
		model_strategy->setup_suffix(model);
		model_strategy->setup_sys_size(model);
		model_strategy->setup_aux_data(model);
		model_strategy->setup_period(model);
		model_strategy->setup_hamiltonian(model);
		model_strategy->setup_hamiltonian_drv(model);
		model_strategy->setup_dissipators(model);
		model_strategy->setup_lindbladian(model);
		model_strategy->setup_lindbladian_drv(model);
		
		model.save_data();
		model.log_setup_info();
	}
};
