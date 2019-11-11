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
		const runge_kutta4_stepper rk4_stepper;
		boost::numeric::odeint::integrate_times(
			rk4_stepper,
			MBLSystem(model),
			start_state,
			times,
			step,
			MBLObserver(model, times, start_state)
		);
	}
};
