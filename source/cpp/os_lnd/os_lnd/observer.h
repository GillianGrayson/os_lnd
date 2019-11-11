#pragma once
#include "model.h"

struct BaseObserver
{
	Model& model;
	std::vector<double>& times;
	
	Eigen::VectorXcd base_state;

	BaseObserver(Model& model, std::vector<double>& times, Eigen::VectorXcd& base_state) : model(model), times(times), base_state(base_state)
	{
	}

	virtual ~BaseObserver() = default;

	virtual void operator()(const Eigen::VectorXcd& x, double t) = 0;
};
