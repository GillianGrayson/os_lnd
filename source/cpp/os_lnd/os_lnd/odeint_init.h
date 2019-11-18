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
	const double start_observed_period = model.ini.GetReal("odeint", "start_observed_period", 0.0);
	const double finish_observed_period = model.ini.GetReal("odeint", "finish_observed_period", 0.0);
	const int num_time_points = model.ini.GetInteger("odeint", "num_time_points", 0);

	std::vector<double> times;
	if (start_observed_period > 0.0)
	{
		times.push_back(start_time);
	}

	if (dump_type == "linear")
	{
		const double shift = (finish_observed_period - start_observed_period) / double(num_time_points - 1);
		for (auto t_id = 0; t_id < num_time_points; t_id++)
		{
			times.push_back(start_time + (start_observed_period + double(t_id) * shift) * model.period);
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
