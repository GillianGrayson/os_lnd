#pragma once
#include "model.h"

struct BaseSystem
{
	Model model;

	BaseSystem(Model& model) : model(model)
	{
	}

	virtual ~BaseSystem() = default;

	virtual void operator()(const Eigen::VectorXcd& x, Eigen::VectorXcd& dxdt, const double t) = 0;
};

struct MBLSystem : BaseSystem
{
	MBLSystem(Model& model) : BaseSystem(model)
	{
	}

	void operator()(const Eigen::VectorXcd& x, Eigen::VectorXcd& dxdt, const double t) override
	{
		dxdt.noalias() = model.lindbladian * x;
	}
};
