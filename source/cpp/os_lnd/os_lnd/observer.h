#pragma once
#include "model.h"
#include <iostream>
#include <iomanip>

struct BaseObserver
{
	Model model;

	BaseObserver(Model& model) : model(model)
	{
	}

	virtual ~BaseObserver() = default;

	virtual void operator()(const Eigen::VectorXcd& x, double t) = 0;
};

struct MBLObserver : BaseObserver
{
	Eigen::VectorXcd prev_state;

	MBLObserver(Model& model, const Eigen::VectorXcd& state) : BaseObserver(model), prev_state(state)
	{
	}

	void operator()(const Eigen::VectorXcd& x, double t) override
	{
		const int precision = model.ini.GetInteger("global", "save_precision", 0);

		Eigen::VectorXcd state_diff = x - prev_state;
		std::cout << "time = " << t << "    diff = " << std::setprecision(precision) << std::scientific << state_diff.norm() << std::endl;
		prev_state = x;
	}
};
