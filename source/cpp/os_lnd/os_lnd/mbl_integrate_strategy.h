#pragma once
#include "integrate_strategy.h"
#include "model.h"
#include <boost/numeric/odeint.hpp>
#include "odeint_init.h"
#include "mbl_system.h"
#include "mbl_observer.h"
#include "save.h"


struct MBLIntegrateStrategy : IntegrateStrategy
{
	MBLIntegrateStrategy(Model& model, std::vector<double>& times, double& step, Eigen::VectorXcd& start_state) : IntegrateStrategy(model, times, step, start_state)
	{
	}

	void integrate_times_rk4() override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);

		const runge_kutta4_stepper rk4_stepper;
		boost::numeric::odeint::integrate_times(
			rk4_stepper,
			MBLSystem(model),
			start_state,
			times,
			step,
			MBLObserver(model, times, start_state)
		);

		model.rho = Eigen::Map<Eigen::MatrixXcd>(start_state.data(), model.sys_size, model.sys_size);

		const auto fn = "rho_mtx" + model.suffix;
		save_dense_mtx(model.rho, fn, save_precision);
	}
};
