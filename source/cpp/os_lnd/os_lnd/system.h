#pragma once
#include "model.h"
#include "init.h"

struct BaseSystem
{
	Model& model;

	BaseSystem(Model& model) : model(model)
	{
	}

	virtual ~BaseSystem() = default;

	virtual void operator()(const Eigen::VectorXcd& x, Eigen::VectorXcd& dxdt, const double t) = 0;
};

