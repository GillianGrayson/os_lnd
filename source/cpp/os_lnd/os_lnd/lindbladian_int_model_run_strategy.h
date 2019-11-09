#pragma once
#include "model_run_strategy.h"
#include <Eigen/SparseLU>
#include "save.h"
#include <iostream>
#include <boost/numeric/odeint.hpp>
#include "odeint_init.h"
#include "system.h"
#include "observer.h"

struct LindbladianIntModelRunStrategy : ModelRunStrategy
{
	void run(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const int start_state_id = model.ini.GetInteger("odeint", "start_state_id", 0);
		const double abs_err = model.ini.GetReal("odeint", "abs_err", 0);
		const double rel_err = model.ini.GetReal("odeint", "rel_err", 0);

		Eigen::VectorXcd start_state = Eigen::VectorXcd::Zero(model.sys_size * model.sys_size);
		start_state[start_state_id * model.sys_size + start_state_id] = std::complex<double>(1.0, 0.0);

		const boost::numeric::odeint::runge_kutta4<Eigen::VectorXcd> runge_kutta4_stepper;
		boost::numeric::odeint::integrate_const(runge_kutta4_stepper, MBLSystem(model), start_state, 0.0, 1000.0, 0.001, MBLObserver(model, start_state));

		//boost::numeric::odeint::integrate_adaptive(
		//	boost::numeric::odeint::make_controlled(abs_err, rel_err, error_stepper_type()),
		//	MBLSystem(model),
		//	start_state,
		//	0.0,
		//	1000.0,
		//	0.0001,
		//	MBLObserver(model, start_state)
		//);

		int a = 0;
	}
};
