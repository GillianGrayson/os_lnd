#pragma once
#include "model.h"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include <Eigen/Core>

typedef boost::numeric::odeint::runge_kutta_dopri5<Eigen::VectorXcd> error_stepper_type;
typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
typedef boost::numeric::odeint::runge_kutta4<Eigen::VectorXcd> runge_kutta4_stepper;


inline std::vector<double> get_times_vector(Model& model)
{
	const std::string dump_type = model.ini.Get("odeint", "dump_type", "unknown");
	const double start_time = model.ini.GetReal("odeint", "start_time", 0.0);
	const double finish_time = model.ini.GetReal("odeint", "finish_time", 0.0);
	const int num_time_points = model.ini.GetInteger("odeint", "num_time_points", 0);

	std::vector<double> times;
	if (start_time > 0.0)
	{
		times.push_back(0.0);
	}

	if (dump_type == "linear")
	{
		const double shift = (finish_time - start_time) / double(num_time_points - 1);
		for (auto t_id = 0; t_id < num_time_points; t_id++)
		{
			times.push_back(start_time + double(t_id) * shift);
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
