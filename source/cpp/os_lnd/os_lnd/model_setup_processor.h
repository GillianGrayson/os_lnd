#pragma once
#include "mbl_model_setup_strategy.h"
#include "routines.h"


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
			model.throw_error("Unsupported model strategy");
		}
	}

	void process_model(Model& model) const
	{
		model_strategy->setup_suffix(model);
		model_strategy->setup_sys_size(model);
		model_strategy->setup_hamiltonian(model);
		model_strategy->setup_hamiltonian_drv(model);
		model_strategy->setup_dissipators(model);
		model_strategy->setup_lindbladian(model);
		model_strategy->setup_lindbladian_drv(model);

		model.log_message(fmt::format("sys_size = {}", model.sys_size));
		model.log_message(fmt::format("Number of non-zero elements in lindbladian = {}", model.lindbladian.nonZeros()));
		model.log_message(fmt::format("Part of non-zero elements in lindbladian = {:.16e}\n", double(model.lindbladian.nonZeros()) / (std::pow(double(model.sys_size), 4.0))));
	}
};