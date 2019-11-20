#pragma once
#include "model.h"
#include <ghc/filesystem.hpp>
#include <stdio.h>


struct BaseObserver
{
	Model& model;
	std::vector<double>& times;
	
	Eigen::VectorXcd base_state;
	std::vector<double> passed_times;

	bool dump_progress;

	BaseObserver(Model& model, std::vector<double>& times, Eigen::VectorXcd& base_state) : model(model), times(times), base_state(base_state)
	{
		dump_progress = model.ini.GetBoolean("odeint", "dump_progress", false);
	}

	virtual ~BaseObserver() = default;

	virtual void operator()(const Eigen::VectorXcd& x, double t) = 0;

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
