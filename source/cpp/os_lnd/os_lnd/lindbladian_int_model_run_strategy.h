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
		const int start_state_id = model.ini.GetInteger("int", "start_state_id", 0);
		const double abs_err = model.ini.GetReal("int", "abs_err", 0);
		const double rel_err = model.ini.GetReal("int", "rel_err", 0);

		Eigen::VectorXcd start_state(std::complex<double>(0.0, 0.0), model.sys_size * model.sys_size);
		start_state[start_state_id].real(1.0);

		boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(abs_err, rel_err), MBLSystem(model), start_state, 0.0, 1000.0, 0.0001, MBLObserver(model, start_state));


		int a = 0;
	}
};
