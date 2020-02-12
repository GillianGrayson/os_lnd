#pragma once
#include "integrate_strategy.h"
#include "model.h"
#include <boost/numeric/odeint.hpp>
#include "odeint_init.h"
#include "mbl_system.h"
#include "mbl_observer.h"


struct MBLIntegrateStrategy : IntegrateStrategy
{
	MBLIntegrateStrategy(Model& model, std::vector<double>& times, double& step, Eigen::VectorXcd& start_state) : IntegrateStrategy(model, times, step, start_state)
	{
	}

	void integrate_times_rk4() override
	{
		MBLSystem system(model);
		
		MBLObserver observer(
			model,
			times,
			start_state
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

		if (!observer.dump_progress)
		{
			double t = observer.passed_times.back();
			
			observer.rewrite_observables("times", observer.passed_times, t, t);
			observer.rewrite_observables("diffs", observer.diffs, t, t);
			observer.dump_current_state(t, observer.diffs.back());

			observer.rewrite_observables("ratios", observer.ratios, t, t);
			observer.rewrite_observables("ees", observer.ees, t, t);
			observer.rewrite_observables("imbalances", observer.imbalances, t, t);
		}
	}
};
