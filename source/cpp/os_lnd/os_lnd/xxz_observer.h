#pragma once
#include "observer.h"
#include "xxz_model_strategy.h"


struct XXZObserver : BaseObserver
{
	XXZModelStrategy model_strategy;
	std::vector<double> znds;
	std::vector<double> vaks;

	XXZObserver(
		Model& model,
		std::vector<double>& times,
		Eigen::VectorXcd& base_state
	) : BaseObserver(model, times, base_state)
	{
		model_strategy = XXZModelStrategy();
		model_strategy.setup_aux_data(model);
	}

	void operator()(const Eigen::VectorXcd& x, double t) override
	{
		process_observables_basic(x, t);

		double znd = model_strategy.get_quantity_znd(model);
		znds.push_back(znd);
		double vak = model_strategy.get_quantity_vak(model);
		vaks.push_back(vak);

		if (dump_progress || is_last_time(t))
		{
			rewrite_observables("znds", znds, t_pre, t);
			rewrite_observables("vaks", vaks, t_pre, t);
		}
	}
};
