#pragma once
#include "mbl_model_strategy.h"
#include "dimer_model_strategy.h"
#include "super_decoh_model_strategy.h"
#include "xxz_model_strategy.h"

struct ModelProcessor
{
	std::unique_ptr<ModelStrategy> model_strategy;

	void set_strategy(Model& model)
	{
		const std::string system = model.ini.Get("global", "system", "unknown");
		const std::string run_type = model.ini.Get("global", "run_type", "unknown");
		if (system == "mbl")
		{
			if (run_type == "regular")
			{
				model_strategy = std::make_unique<MBLModelStrategy>();
			}
			else
			{
				model.throw_error(fmt::format("Unsupported run_type for {s}", system));
			}
		}
		else if (system == "xxz")
		{
			if (run_type == "regular" || run_type == "serial")
			{
				model_strategy = std::make_unique<XXZModelStrategy>();
			}
			else
			{
				model.throw_error(fmt::format("Unsupported run_type for {s}", system));
			}
		}
		else if (system == "dimer")
		{
			if (run_type == "regular")
			{
				model_strategy = std::make_unique<DimerModelStrategy>();
			}
			else
			{
				model.throw_error(fmt::format("Unsupported run_type for {s}", system));
			}
		}
		else if (system == "super_decoh")
		{
			if (run_type == "regular" || run_type == "serial")
			{
				model_strategy = std::make_unique<SuperDecohModelStrategy>();
			}
			else
			{
				model.throw_error(fmt::format("Unsupported run_type for {s}", system));
			}
		}
		else
		{
			model.throw_error("Unsupported model_strategy");
		}
	}

	void set_serial_features(
		Model& model,
		std::map<std::string, std::vector<double>>& features_double, 
		std::map<std::string, std::vector<std::complex<double>>>& features_complex)
	{
		model_strategy->setup_serial_data(model, features_double, features_complex);
	}

	void fill_serial_features(
		Model& model,
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex)
	{
		model_strategy->fill_serial_data(model, features_double, features_complex);
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
		model_strategy->setup_lindbladians_drv(model);
		
		model.save_data();
		model.log_setup_info();
	}
};
