#pragma once
#include "model.h"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include <Eigen/Core>

typedef boost::numeric::odeint::runge_kutta_dopri5<Eigen::VectorXcd> error_stepper_type;
typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
typedef boost::numeric::odeint::runge_kutta4<Eigen::VectorXcd> runge_kutta4_stepper;


inline std::vector<double> get_times_vector(Model& model, double start_time)
{
	const std::string dump_type = model.ini.Get("odeint", "dump_type", "unknown");
	const double total_num_periods = model.ini.GetReal("odeint", "total_num_periods", 0.0);
	const double current_num_periods = model.ini.GetReal("odeint", "current_num_periods", 0.0);
	const int current_num_time_points = model.ini.GetInteger("odeint", "current_num_time_points", 0);

	std::vector<double> times;

	if (dump_type == "linear")
	{
		const double shift = current_num_periods / double(current_num_time_points - 1);
		
		if (start_time > (total_num_periods - shift * 0.5) * model.period)
		{
			model.throw_error("Integration already done");
		}
		
		for (auto t_id = 0; t_id < current_num_time_points; t_id++)
		{
			double time = start_time + (double(t_id) * shift) * model.period;
			if (time < (total_num_periods + shift * 0.5) * model.period)
			{
				times.push_back(time);
			}
			else
			{
				break;
			}
		}
	}
	else if (dump_type == "log")
	{
		model.throw_error("Unsupported dump_type");
	}
	else
	{
		model.throw_error("Unsupported dump_type");
	}

	return times;
}
