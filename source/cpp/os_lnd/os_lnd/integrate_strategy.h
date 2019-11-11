#pragma once
#include "model.h"
#include "system.h"

struct IntegrateStrategy
{
	Model model;
	std::vector<double> times;
	double step;
	Eigen::VectorXcd start_state;

	IntegrateStrategy(Model& model, std::vector<double>& times, double& step, Eigen::VectorXcd& start_state) : model(model), times(times), step(step), start_state(start_state)
	{
	}

	virtual ~IntegrateStrategy() = default;

	virtual void integrate_times_rk4() = 0;
};
