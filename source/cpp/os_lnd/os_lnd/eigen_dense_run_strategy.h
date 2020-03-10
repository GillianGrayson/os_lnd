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

		model.log_message("Lindbladians eigen...");
		
		Eigen::ComplexEigenSolver<ds_mtx> es;
		es.compute(model.lindbladian_dense, true);

		auto tmp = es.eigenvalues();
		std::vector<std::complex<double>> eigen_values(tmp.data(), tmp.data() + tmp.rows() * tmp.cols());

		auto fn = "lindbladian_evals" + model.suffix;
		save_vector(eigen_values, fn, save_precision);
		
		/*auto min_eval_index = std::min_element(eigen_values.begin(), eigen_values.end()) - eigen_values.begin();

		model.log_message(fmt::format("min_eval_index = {:d}", min_eval_index));
		model.log_message(fmt::format("min_eval_index = {:16e} + {:16e} i", eigen_values[min_eval_index].real(), eigen_values[min_eval_index].imag()));

		auto rho_vec = es.eigenvectors().col(min_eval_index);*/

		/*model.rho = Eigen::Map<Eigen::MatrixXcd>(rho_vec.data(), model.sys_size, model.sys_size);
		std::complex<double> trace_rho = model.rho.trace();
		model.log_message(fmt::format("trace_rho = {:16e} + {:16e} i", trace_rho.real(), trace_rho.imag()));
		model.rho = model.rho / trace_rho;
		trace_rho = model.rho.trace();
		model.log_message(fmt::format("trace_rho = {:16e} + {:16e} i", trace_rho.real(), trace_rho.imag()));
		
		fn = "rho_mtx" + model.suffix;
		save_dense_mtx(model.rho, fn, save_precision);*/
	}
};