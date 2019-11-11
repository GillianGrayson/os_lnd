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
	Model model;
	Eigen::VectorXcd prev_state;

	MBLObserver(Model& model, const Eigen::VectorXcd& state) : BaseObserver(model), model(model), prev_state(state)
	{
	}

	void operator()(const Eigen::VectorXcd& x, double t) override
	{
		const int precision = model.ini.GetInteger("global", "save_precision", 0);

		const Eigen::VectorXcd state_diff = x - prev_state;
		prev_state = x;

		model.log_time_duration();
		model.log_message(fmt::format("time = {:.16e}", t));
		model.log_message(fmt::format("diff = {:.16e}\n", state_diff.norm()));
	}
};
