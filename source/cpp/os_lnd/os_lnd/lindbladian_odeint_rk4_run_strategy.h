#pragma once
#include "run_strategy.h"
#include "save.h"
#include <boost/numeric/odeint.hpp>
#include "odeint_init.h"
#include "system.h"
#include "observer.h"

struct LindbladianODEIntRK4RunStrategy : RunStrategy
{
	void run(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const int start_state_id = model.ini.GetInteger("odeint", "start_state_id", 0);
		const double step = model.ini.GetReal("odeint", "step", 0.0);
	
		Eigen::VectorXcd start_state = Eigen::VectorXcd::Zero(model.sys_size * model.sys_size);
		start_state[start_state_id * model.sys_size + start_state_id] = std::complex<double>(1.0, 0.0);

		const runge_kutta4_stepper rk4_stepper;
		std::vector<double> times = get_times_vector(model);

		boost::numeric::odeint::integrate_times(
			rk4_stepper,
			MBLSystem(model),
			start_state,
			times,
			step,
			MBLObserver(model, start_state)
		);

		model.rho = Eigen::Map<Eigen::MatrixXcd>(start_state.data(), model.sys_size, model.sys_size);

		const auto fn = "rho_mtx" + model.suffix;
		save_dense_mtx(model.rho, fn, save_precision);
	}
};
