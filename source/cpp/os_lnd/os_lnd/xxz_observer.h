#pragma once
#include "observer.h"
#include "xxz_model_strategy.h"


struct XXZObserver : BaseObserver
{
	XXZModelStrategy model_strategy;
	std::vector<std::vector<double>> js;
	int num_spins;

	XXZObserver(
		Model& model,
		std::vector<double>& times,
		Eigen::VectorXcd& base_state
	) : BaseObserver(model, times, base_state)
	{
		model_strategy = XXZModelStrategy();
		model_strategy.setup_aux_data(model);

		num_spins = model.ini.GetInteger("xxz", "num_spins", 0);
		for (auto spin_id = 0; spin_id < num_spins - 1; spin_id++)
		{
			js.push_back({});
		}
	}

	void operator()(const Eigen::VectorXcd& x, double t) override
	{
		process_observables_basic(x, t);

		std::vector<double> quantities = model_strategy.get_quantities(model);
		for (auto spin_id = 0; spin_id != quantities.size(); spin_id++) 
		{
			js[spin_id].push_back(quantities[spin_id]);
		}

		if (is_dump_now(t))
		{
			for (auto spin_id = 0; spin_id != quantities.size(); spin_id++)
			{
				rewrite_observables("j_" + std::to_string(spin_id), js[spin_id], t_pre, t);
			}
		}
	}

	void fill_serial_features(
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex) override
	{
		for (auto spin_id = 0; spin_id < num_spins - 1; spin_id++)
		{
			features_double["j_" + std::to_string(spin_id)].insert(features_double["j_" + std::to_string(spin_id)].end(), js[spin_id].begin(), js[spin_id].end());
		}
	}
};
