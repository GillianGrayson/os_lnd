#pragma once
#include "model.h"
#include "init.h"

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
		//std::cout << x.size() << std::endl;
		//std::cout << dxdt.size() << std::endl;
		dxdt.noalias() = model.lindbladian * x;
		//std::cout << dxdt.size() << std::endl;
	}
};
