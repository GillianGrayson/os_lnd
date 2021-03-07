#pragma once
#include "observer.h"


struct MBLObserver : BaseObserver
{
	MBLModelStrategy model_strategy;
	std::vector<double> ratios;
	std::vector<double> ees;
	std::vector<double> imbalances;
	
	MBLObserver(
		Model& model,
		std::vector<double>& times,
		Eigen::VectorXcd& base_state
	) : BaseObserver(model, times, base_state)
	{
		model_strategy = MBLModelStrategy();
		model_strategy.setup_aux_data(model);
	}

	void operator()(const Eigen::VectorXcd& x, double t) override
	{		
		process_observables_basic(x, t);
		
		double ratio = model_strategy.get_ratio(model);
		ratios.push_back(ratio);

		double ee = model_strategy.get_entanglement_entropy(model);
		ees.push_back(ee);

		double imbalance = model_strategy.get_imbalance(model);
		imbalances.push_back(imbalance);
		
		if (dump_progress || is_last_time(t))
		{
			rewrite_observables("ratios", ratios, t_pre, t);
			rewrite_observables("ees", ees, t_pre, t);
			rewrite_observables("imbalances", imbalances, t_pre, t);
		}
	}

	void fill_serial_features(
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex) override
	{
	}
};
