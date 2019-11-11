#pragma once
#include "mbl_setup_strategy.h"
#include "routines.h"


struct SetupProcessor
{
	std::unique_ptr<SetupStrategy> setup_strategy;

	void set_model_strategy(Model& model)
	{
		const std::string system = model.ini.Get("global", "system", "unknown");
		if (system == "mbl")
		{
			setup_strategy = std::make_unique<MBLSetupStrategy>();
		}
		else
		{
			model.throw_error("Unsupported setup_strategy");
		}
	}

	void process_model(Model& model) const
	{
		setup_strategy->setup_suffix(model);
		setup_strategy->setup_sys_size(model);
		setup_strategy->setup_hamiltonian(model);
		setup_strategy->setup_hamiltonian_drv(model);
		setup_strategy->setup_dissipators(model);
		setup_strategy->setup_lindbladian(model);
		setup_strategy->setup_lindbladian_drv(model);

		model.log_setup_info();
	}
};