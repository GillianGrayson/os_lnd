#pragma once
#include "mbl_integrate_strategy.h"
#include "dimer_integrate_strategy.h"
#include "xxz_integrate_strategy.h"


struct IntegrateProcessor
{
	std::unique_ptr<IntegrateStrategy> integrate_strategy;

	void set_strategy(Model& model, std::vector<double>& times, double& step, Eigen::VectorXcd& start_state)
	{
		const std::string system = model.ini.Get("global", "system", "unknown");
		if (system == "mbl")
		{
			integrate_strategy = std::make_unique<MBLIntegrateStrategy>(model, times, step, start_state);
		}
		else if (system == "dimer")
		{
			integrate_strategy = std::make_unique<DimerIntegrateStrategy>(model, times, step, start_state);
		}
		else if (system == "xxz")
		{
			integrate_strategy = std::make_unique<XXZIntegrateStrategy>(model, times, step, start_state);
		}
		else
		{
			model.throw_error("Unsupported system");
		}
	}

	void process() const
	{
		integrate_strategy->integrate_times_rk4();
	}
};
