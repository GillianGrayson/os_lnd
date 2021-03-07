#pragma once
#include "run_strategy.h"
#include <Eigen/Eigenvalues>
#include "save.h"
#include <algorithm>

struct EigenDenseRunStrategy : RunStrategy
{
	void eigen(Model& model)
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const auto save_lindbladian_evals = model.ini.GetBoolean("global", "save_lindbladian_evals", false);
		const auto save_rho = model.ini.GetBoolean("global", "save_rho", false);
		const auto save_rho_evals = model.ini.GetBoolean("global", "save_rho_evals", false);

		auto is_sp_empty = (model.lindbladian.outerSize() > 0) ? false : true;
		auto is_ds_empty = (model.lindbladian_dense.outerSize() > 0) ? false : true;

		if (is_sp_empty && is_ds_empty)
		{
			model.throw_error("Lindbladian is empty");
		}

		if (is_ds_empty && !is_sp_empty)
		{
			model.lindbladian_dense = ds_mtx(model.lindbladian);
		}

		model.log_time_duration();
		model.log_message("Lindbladians eigen...");

		Eigen::ComplexEigenSolver<ds_mtx> es;
		es.compute(model.lindbladian_dense, true);

		auto lind_evals_tmp = es.eigenvalues();
		model.lindbladian_evals = std::vector<std::complex<double>>(lind_evals_tmp.data(), lind_evals_tmp.data() + lind_evals_tmp.rows() * lind_evals_tmp.cols());
		if (save_lindbladian_evals)
		{
			auto fn = "lindbladian_evals" + model.suffix;
			save_vector(model.lindbladian_evals, fn, save_precision);
		}

		Eigen::VectorXcd zeros = Eigen::VectorXcd::Zero(model.sys_size);
		std::vector<double> evec_sub_diag_norms(model.sys_size * model.sys_size);
		for (auto i = 0; i < model.sys_size * model.sys_size; i++)
		{
			Eigen::MatrixXcd evec = es.eigenvectors().col(i);
			Eigen::MatrixXcd evec_mtx = Eigen::Map<Eigen::MatrixXcd>(evec.data(), model.sys_size, model.sys_size);
			evec_mtx.diagonal() = zeros;
			evec_sub_diag_norms[i] = evec_mtx.norm();
		}
		if (save_lindbladian_evals)
		{
			auto fn = "evec_sub_diag_norms" + model.suffix;
			save_vector(evec_sub_diag_norms, fn, save_precision);
		}

		std::vector<double> abs_evals(model.sys_size * model.sys_size);
		for (auto i = 0; i < model.sys_size * model.sys_size; i++)
		{
			abs_evals[i] = std::abs(model.lindbladian_evals[i]);
		}

		auto min_eval_index = std::distance(abs_evals.begin(), std::min_element(abs_evals.begin(), abs_evals.end()));

		model.log_message(fmt::format("min_eval_index = {:d}", min_eval_index));
		model.log_message(fmt::format("min_eval = {:16e} + {:16e} i", model.lindbladian_evals[min_eval_index].real(), model.lindbladian_evals[min_eval_index].imag()));

		Eigen::MatrixXcd rho_vec = es.eigenvectors().col(min_eval_index);

		model.rho = Eigen::Map<Eigen::MatrixXcd>(rho_vec.data(), model.sys_size, model.sys_size);
		std::complex<double> trace_rho = model.rho.trace();
		model.log_message(fmt::format("trace_rho = {:16e} + {:16e} i", trace_rho.real(), trace_rho.imag()));
		model.rho = model.rho / trace_rho;
		trace_rho = model.rho.trace();
		model.log_message(fmt::format("trace_rho = {:16e} + {:16e} i", trace_rho.real(), trace_rho.imag()));
		model.log_time_duration();
		if (save_rho)
		{
			auto fn = "rho_mtx" + model.suffix;
			save_dense_mtx(model.rho, fn, save_precision);
		}
		es.compute(model.rho, true);
		auto rho_evals_tmp = es.eigenvalues();
		model.rho_evals = std::vector<std::complex<double>>(rho_evals_tmp.data(), rho_evals_tmp.data() + rho_evals_tmp.rows() * rho_evals_tmp.cols());
		if (save_rho_evals)
		{
			auto fn = "rho_evals" + model.suffix;
			save_vector(model.rho_evals, fn, save_precision);
		}
	}
	
	void run(Model& model) override
	{
		eigen(model);
		
		ModelProcessor model_processor;
		model_processor.set_strategy(model);
		model_processor.init_model(model);
		model_processor.release_observables(model);
		model.log_time_duration();
	}

	void run_serial(
		Model& model,
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex) override
	{
		eigen(model);
		
		ModelProcessor model_processor;
		model_processor.set_strategy(model);
		model_processor.init_model(model);
		model_processor.fill_serial_features(model, features_double, features_complex);
		model.log_time_duration();
	}
};