#pragma once
#include "system.h"

struct PImplSystem
{
	std::shared_ptr<BaseSystem> system;

	PImplSystem()
	{
	}

	PImplSystem(std::shared_ptr<BaseSystem> sys) : system(std::move(sys))
	{
	}

	void operator()(const Eigen::VectorXcd& x, Eigen::VectorXcd& dxdt, const double t)
	{
		system->operator()(x, dxdt, t);
	}
};
