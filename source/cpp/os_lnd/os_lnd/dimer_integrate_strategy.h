#pragma once
#include "integrate_strategy.h"
#include "model.h"
#include <boost/numeric/odeint.hpp>
#include "odeint_init.h"
#include "mbl_system.h"
#include "mbl_observer.h"
#include "save.h"


struct DimerIntegrateStrategy : IntegrateStrategy
{
	DimerIntegrateStrategy(Model& model, std::vector<double>& times, double& step, Eigen::VectorXcd& start_state) : IntegrateStrategy(model, times, step, start_state)
	{
	}

	void integrate_times_rk4() override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);

		MBLSystem system(model);

		std::vector<double> diffs;
		MBLObserver observer(
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

		model.rho = Eigen::Map<Eigen::MatrixXcd>(start_state.data(), model.sys_size, model.sys_size);

		auto fn = "rho_mtx" + model.suffix;
		save_dense_mtx(model.rho, fn, save_precision);

		fn = "diffs" + model.suffix;
		save_vector(diffs, fn, save_precision);

		fn = "times" + model.suffix;
		save_vector(times, fn, save_precision);
	}
};
