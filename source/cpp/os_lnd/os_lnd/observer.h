#pragma once
#include "model.h"
#include <ghc/filesystem.hpp>
#include <stdio.h>


struct BaseObserver
{
	Model& model;
	std::vector<double>& times;
	Eigen::VectorXcd base_state; // not a reference
	
	bool dump_progress;
	std::vector<double> passed_times;
	std::vector<double> diffs;

	BaseObserver(Model& model, std::vector<double>& times, Eigen::VectorXcd& base_state) : model(model), times(times), base_state(base_state)
	{
		dump_progress = model.ini.GetBoolean("odeint", "dump_progress", false);
	}

	virtual ~BaseObserver() = default;

	virtual void operator()(const Eigen::VectorXcd& x, double t) = 0;

	void process_observables_basic(const Eigen::VectorXcd& x, double t)
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

	std::string get_suffix(const double t) const
	{
		std::stringstream fns;
		fns << "_times(" << std::setprecision(2) << std::scientific << times[0] << "_" << t << ")";
		return fns.str();
	}

	template <typename T>
	void rewrite_observables(const std::string& name, const std::vector<T>& observable, const double t_pre, const double t_now)
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		
		std::string fn;

		fn = name + get_suffix(t_pre) + model.suffix;
		if (ghc::filesystem::exists(fn))
		{
			remove(fn.c_str());
		}

		fn = name + get_suffix(t_now) + model.suffix;
		save_vector(observable, fn, save_precision);
	}

	void dump_current_state(const double t, const double diff)
	{
		if (dump_progress)
		{
			const int save_precision = model.ini.GetInteger("global", "save_precision", 0);

			model.rho = Eigen::Map<Eigen::MatrixXcd>(base_state.data(), model.sys_size, model.sys_size);
			auto fn = "rho_mtx" + model.suffix;
			save_dense_mtx(model.rho, fn, save_precision);

			std::vector<double> curr_dump;
			curr_dump.push_back(t);
			curr_dump.push_back(diff);

			fn = "curr_dump" + model.suffix;
			save_vector(curr_dump, fn, save_precision);
		}
	}
};
