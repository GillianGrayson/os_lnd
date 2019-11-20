#pragma once
#include "integrate_strategy.h"
#include "model.h"
#include <boost/numeric/odeint.hpp>
#include "odeint_init.h"
#include "dimer_system.h"
#include "dimer_observer.h"


struct DimerIntegrateStrategy : IntegrateStrategy
{
	DimerIntegrateStrategy(Model& model, std::vector<double>& times, double& step, Eigen::VectorXcd& start_state) : IntegrateStrategy(model, times, step, start_state)
	{
	}

	void integrate_times_rk4() override
	{
		DimerSystem system(model);

		std::vector<double> diffs;
		DimerObserver observer(
			model,
			times,
			start_state,
			diffs
		);

		const runge_kutta4_stepper rk4_stepper;
		boost::numeric::odeint::integrate_times(
			rk4_stepper,
			system,
			start_state,
			times,
			step,
			observer
		);
	}
};
