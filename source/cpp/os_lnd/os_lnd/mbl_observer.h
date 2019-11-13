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
		if (t > std::numeric_limits<double>::epsilon())
		{
			model.log_time_duration();
			
			const Eigen::VectorXcd state_diff = x - base_state;
			base_state = x;

			double diff = state_diff.norm();

			if (std::isnan(diff))
			{
				model.throw_error("Integration failed. Try to decrease integration step");
			}

			diffs.push_back(diff);

			model.log_message(fmt::format("time = {:.16e}", t));
			model.log_message(fmt::format("diff = {:.16e}\n", diff));

			dump_current_state(t, diff);
		}
	}
};
