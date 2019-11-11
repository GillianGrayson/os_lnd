#pragma once
#include "model.h"
#include "save.h"
#include <boost/numeric/odeint.hpp>
#include "odeint_init.h"
#include "system.h"
#include "observer.h"

struct IntegrateStrategy
{
	Model model;
	std::vector<double> times;
	Eigen::VectorXcd start_state;

	IntegrateStrategy(Model& model, std::vector<double>& times, Eigen::VectorXcd& base_state) : model(model), times(times), start_state(base_state)
	{
	}

	virtual ~IntegrateStrategy() = default;

	virtual void integrate_times_rk4() = 0;
};

struct MBLIntegrator : IntegrateStrategy
{
	MBLIntegrator(Model& model, std::vector<double>& times, Eigen::VectorXcd& base_state) : IntegrateStrategy(model, times, base_state)
	{
	}
		
	void integrate_times_rk4() override
	{
		const double step = model.ini.GetReal("odeint", "step", 0.0);

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
	}
};
