#pragma once
#include "observer.h"


struct MBLObserver : BaseObserver
{	
	MBLObserver(
		Model& model,
		std::vector<double>& times,
		Eigen::VectorXcd& base_state
	) : BaseObserver(model, times, base_state)
	{
	}

	void operator()(const Eigen::VectorXcd& x, double t) override
	{
		process_observables_basic(x, t);
	}
};
