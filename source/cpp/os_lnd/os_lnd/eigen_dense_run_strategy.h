#pragma once
#include "run_strategy.h"
#include <Eigen/Eigenvalues>
#include "save.h"
#include <algorithm>

struct EigenDenseRunStrategy : RunStrategy
{
	void run(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const auto save_rho = model.ini.GetBoolean("global", "save_rho", false);

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
		std::vector<std::complex<double>> lind_evals(lind_evals_tmp.data(), lind_evals_tmp.data() + lind_evals_tmp.rows() * lind_evals_tmp.cols());
		auto fn = "lindbladian_evals" + model.suffix;
		save_vector(lind_evals, fn, save_precision);

		Eigen::VectorXcd zeros = Eigen::VectorXcd::Zero(model.sys_size);
		std::vector<double> evec_sub_diag_norms(model.sys_size * model.sys_size);
		for (auto i = 0; i < model.sys_size * model.sys_size; i++)
		{
			Eigen::MatrixXcd evec = es.eigenvectors().col(i);
			Eigen::MatrixXcd evec_mtx = Eigen::Map<Eigen::MatrixXcd>(evec.data(), model.sys_size, model.sys_size);
			evec_mtx.diagonal() = zeros;
			evec_sub_diag_norms[i] = evec_mtx.norm();
		}
		fn = "evec_sub_diag_norms" + model.suffix;
		save_vector(evec_sub_diag_norms, fn, save_precision);

		std::vector<double> abs_evals(model.sys_size * model.sys_size);
		for (auto i = 0; i < model.sys_size * model.sys_size; i++)
		{
			abs_evals[i] = std::abs(lind_evals[i]);
		}

		auto min_eval_index = std::distance(abs_evals.begin(), std::min_element(abs_evals.begin(), abs_evals.end()));

		model.log_message(fmt::format("min_eval_index = {:d}", min_eval_index));
		model.log_message(fmt::format("min_eval = {:16e} + {:16e} i", lind_evals[min_eval_index].real(), lind_evals[min_eval_index].imag()));

		Eigen::MatrixXcd rho_vec = es.eigenvectors().col(min_eval_index);

		model.rho = Eigen::Map<Eigen::MatrixXcd>(rho_vec.data(), model.sys_size, model.sys_size);
		std::complex<double> trace_rho = model.rho.trace();
		model.log_message(fmt::format("trace_rho = {:16e} + {:16e} i", trace_rho.real(), trace_rho.imag()));
		model.rho = model.rho / trace_rho;
		trace_rho = model.rho.trace();
		model.log_message(fmt::format("trace_rho = {:16e} + {:16e} i", trace_rho.real(), trace_rho.imag()));
		model.log_time_duration();
		fn = "rho_mtx" + model.suffix;
		if (save_rho)
		{
			save_dense_mtx(model.rho, fn, save_precision);
		}
		es.compute(model.rho, true);
		auto rho_evals_tmp = es.eigenvalues();
		std::vector<std::complex<double>> rho_evals(rho_evals_tmp.data(), rho_evals_tmp.data() + rho_evals_tmp.rows() * rho_evals_tmp.cols());
		fn = "rho_evals" + model.suffix;
		save_vector(rho_evals, fn, save_precision);

		ModelProcessor model_processor;
		model_processor.set_strategy(model);
		model_processor.init_model(model);
		model_processor.release_observables(model);
		model.log_time_duration();
	}
};