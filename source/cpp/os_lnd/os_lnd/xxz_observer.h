#pragma once
#include "observer.h"
#include "xxz_model_strategy.h"


struct XXZObserver : BaseObserver
{
	XXZModelStrategy model_strategy;
	std::vector<double> quantities;

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

		double quantity = model_strategy.get_quantity(model);
		quantities.push_back(quantity);

		if (dump_progress || is_last_time(t))
		{
			rewrite_observables("quantities", quantities, t_pre, t);
		}
	}
};
