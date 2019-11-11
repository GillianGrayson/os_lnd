#pragma once
#include "run_strategy.h"
#include "save.h"
#include <boost/numeric/odeint.hpp>
#include "odeint_init.h"
#include "integrate_processor.h"

struct LindbladianODEIntRK4RunStrategy : RunStrategy
{
	void run(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const int start_state_id = model.ini.GetInteger("odeint", "start_state_id", 0);
		double step = model.ini.GetReal("odeint", "step", 0.0);
		
		std::vector<double> times = get_times_vector(model);

		Eigen::VectorXcd start_state = Eigen::VectorXcd::Zero(model.sys_size * model.sys_size);
		start_state[start_state_id * model.sys_size + start_state_id] = std::complex<double>(1.0, 0.0);

		IntegrateProcessor integrate_processor;
		integrate_processor.set_strategy(model, times, step, start_state);
		integrate_processor.process();
	}
};
