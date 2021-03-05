#pragma once
#include "run_strategy.h"
#include "load.h"
#include <boost/numeric/odeint.hpp>
#include "odeint_init.h"
#include "integrate_processor.h"
#include <ghc/filesystem.hpp>

struct ODEIntRK4RunStrategy : RunStrategy
{
	void run(Model& model) override
	{
		const int start_state_id = model.ini.GetInteger("odeint", "start_state_id", 0);
		bool is_continue = model.ini.GetBoolean("odeint", "continue", false);
		std::string continue_path = model.ini.Get("odeint", "continue_path", "");
		double step = model.ini.GetReal("odeint", "step", 0.0);

		std::vector<double> times;
		Eigen::VectorXcd start_state;
		auto fn_curr_dump = continue_path + "curr_dump" + model.suffix;
		auto fn_rho_mtx = continue_path + "rho_mtx" + model.suffix;
		
		if (is_continue && ghc::filesystem::exists(fn_curr_dump) && ghc::filesystem::exists(fn_rho_mtx))
		{
			auto fn = continue_path + "curr_dump" + model.suffix;
			std::vector<double> last_dump_info;
			load_vector(last_dump_info, fn);

			fn = continue_path + "rho_mtx" + model.suffix;
			std::vector<std::complex<double>> last_state;
			load_vector(last_state, fn);

			double start_time = last_dump_info[0];
			times = get_times_vector(model, start_time, is_continue);

			start_state = Eigen::Map<Eigen::VectorXcd, Eigen::Unaligned>(last_state.data(), last_state.size());
		}
		else
		{
			times = get_times_vector(model, 0.0, is_continue);

			start_state = Eigen::VectorXcd::Zero(model.sys_size * model.sys_size);
			start_state[start_state_id * model.sys_size + start_state_id] = std::complex<double>(1.0, 0.0);
		}
		
		IntegrateProcessor integrate_processor;
		integrate_processor.set_strategy(model, times, step, start_state);
		integrate_processor.process();
	}
};
