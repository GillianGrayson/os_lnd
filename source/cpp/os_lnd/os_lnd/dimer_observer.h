#pragma once
#include "observer.h"


struct DimerObserver : BaseObserver
{
	std::vector<double>& diffs;
	
	DimerObserver(
		Model& model,
		std::vector<double>& times,
		Eigen::VectorXcd& base_state,
		std::vector<double>& diffs
	) : BaseObserver(model, times, base_state), diffs(diffs)
	{
	}

	void operator()(const Eigen::VectorXcd& x, double t) override
	{
		if (t > times[0] + std::numeric_limits<double>::epsilon())
		{
			model.log_time_duration();
			
			const Eigen::VectorXcd state_diff = x - base_state;
			base_state = x;

			double diff = state_diff.norm();

			if (std::isnan(diff))
			{
				model.throw_error("Integration failed. Try to decrease integration step");
			}

			double t_pre;
			if (!passed_times.empty())
			{
				t_pre = passed_times.back();
			}
			else
			{
				t_pre = 0.0;
			}

			diffs.push_back(diff);
			rewrite_observables("diffs", diffs, t_pre, t);

			passed_times.push_back(t);
			rewrite_observables("times", passed_times, t_pre, t);
			
			model.log_message(fmt::format("time = {:.16e}", t));
			model.log_message(fmt::format("diff = {:.16e}\n", diff));

			dump_current_state(t, diff);
		}
	}
};
