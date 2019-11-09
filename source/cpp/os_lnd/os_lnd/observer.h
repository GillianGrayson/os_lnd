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

		const Eigen::VectorXcd state_diff = x - prev_state;
		prev_state = x;

		const auto run_time = std::chrono::high_resolution_clock::now();
		const auto duration = std::chrono::duration<double>(run_time - model.run_time).count();
		std::cout << std::setprecision(precision) << std::scientific << "time = " << t << "    diff = " << state_diff.norm() << std::endl;
		std::cout << "run_time = " << duration << " seconds" << std::endl << std::endl;
	}
};
