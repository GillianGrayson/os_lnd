#pragma once
#include "observer.h"


struct DimerObserver : BaseObserver
{
	DimerObserver(
		Model& model,
		std::vector<double>& times,
		Eigen::VectorXcd& base_state
	) : BaseObserver(model, times, base_state)
	{
	}

	void operator()(const Eigen::VectorXcd& x, double t) override
	{
		if (t > times[0] + std::numeric_limits<double>::epsilon())
		{
			process_observables_basic(x, t);
		}
	}
};
