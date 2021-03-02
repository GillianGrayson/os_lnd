#pragma once
#include "system.h"


struct XXZSystem : BaseSystem
{
	XXZSystem(Model& model) : BaseSystem(model)
	{
	}

	void operator()(const Eigen::VectorXcd& x, Eigen::VectorXcd& dxdt, const double t) override
	{
		dxdt.noalias() = model.lindbladian * x;
	}
};
