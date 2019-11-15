#pragma once
#include "model.h"
#include "system.h"

struct IntegrateStrategy
{
	Model& model;
	std::vector<double>& times;
	double& step;
	Eigen::VectorXcd& start_state;
	std::string suffix;

	IntegrateStrategy(Model& model, std::vector<double>& times, double& step, Eigen::VectorXcd& start_state) : model(model), times(times), step(step), start_state(start_state)
	{
		std::stringstream fns;
		fns << "_times(" << std::setprecision(2) << std::scientific << times[0] << "_" << times.back() << ")_";
		suffix = fns.str();
	}

	virtual ~IntegrateStrategy() = default;

	virtual void integrate_times_rk4() = 0;
};
