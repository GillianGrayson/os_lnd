#pragma once
#include "observer.h"


struct MBLObserver : BaseObserver
{
	std::vector<double>& diffs;
	
	MBLObserver(
		Model& model,
		std::vector<double>& times,
		Eigen::VectorXcd& base_state,
		std::vector<double>& diffs
	) : BaseObserver(model, times, base_state), diffs(diffs)
	{
	}

	void operator()(const Eigen::VectorXcd& x, double t) override
	{
		const int precision = model.ini.GetInteger("global", "save_precision", 0);

		const Eigen::VectorXcd state_diff = x - base_state;
		base_state = x;

		double diff = state_diff.norm();
		diffs.push_back(diff);

		model.log_time_duration();
		model.log_message(fmt::format("time = {:.16e}", t));
		model.log_message(fmt::format("diff = {:.16e}\n", diff));
	}
};
