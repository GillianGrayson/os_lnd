#pragma once
#include "run_strategy.h"
#include "save.h"
#include <armadillo>


struct LindbladianSmallerEigenVectorRunStrategy : RunStrategy
{
	void run(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);

		std::vector<long long unsigned int> row_ind_vect((model.lindbladian).innerIndexPtr(), (model.lindbladian).innerIndexPtr() + (model.lindbladian).nonZeros());
		std::vector<long long unsigned int> col_ptr_vect((model.lindbladian).outerIndexPtr(), (model.lindbladian).outerIndexPtr() + (model.lindbladian).outerSize() + 1);
		std::vector<std::complex<double>> values_vect((model.lindbladian).valuePtr(), (model.lindbladian).valuePtr() + (model.lindbladian).nonZeros());

		arma::cx_dvec values(values_vect.data(), values_vect.size(), false);
		arma::uvec row_ind(row_ind_vect.data(), row_ind_vect.size(), false);
		arma::uvec col_ptr(col_ptr_vect.data(), col_ptr_vect.size(), false);

		arma::sp_cx_dmat arma_lindbladian(row_ind, col_ptr, values, model.sys_size * model.sys_size, model.sys_size * model.sys_size);
		model.log_message("armadillo lindbladian created");

		arma::cx_vec eigen_values;
		arma::cx_mat eigen_vectors;

		bool result = arma::eigs_gen(eigen_values, eigen_vectors, arma_lindbladian, 1, "sm");

		if (!result)
		{
			model.throw_error("Solving failed");
		}
		else
		{
			model.log_message("Solving complete!");
		}

		//model.rho = Eigen::Map<Eigen::MatrixXcd>(eigen_vectors, model.sys_size, model.sys_size);
		//auto fn = "rho_mtx" + model.suffix;
		//save_dense_mtx(model.rho, fn, save_precision);
	}
};