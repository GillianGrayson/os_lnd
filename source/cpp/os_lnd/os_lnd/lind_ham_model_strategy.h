#pragma once
#include "init.h"
#include "model.h"
#include <gcem.hpp>
#include <mkl.h>
#include <iomanip>
#include "save.h"
#include "model_strategy.h"
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/Eigenvalues>
#include "routines.h"
#include <numeric>

struct LindHamModelStrategy : ModelStrategy
{
	void setup_aux_data(Model& model) override
	{
		model.init_f_basis();
	}

	void setup_suffix(Model& model) override
	{
		const std::string run_type = model.ini.Get("global", "run_type", "unknown");

		int seed;
		if (run_type == "serial")
		{
			seed = std::round(model.serial_state);
		}
		else
		{
			seed = model.ini.GetInteger("lind_ham", "seed", 0);
		}

		const int name_precision = model.ini.GetInteger("global", "name_precision", 0);

		const int N = model.ini.GetInteger("lind_ham", "N", 0);
		const auto alpha = model.ini.GetReal("lind_ham", "alpha", 0.0);

		std::stringstream fns;
		fns << "_N(" << N << ")";
		fns << "_alpha(" << std::setprecision(name_precision) << std::fixed << alpha << ")";

		std::stringstream serial;
		serial << fns.rdbuf();
		if (run_type == "regular")
		{
			fns << "_seed(" << seed << ")";
		}

		fns << ".txt";
		serial << ".txt";

		model.suffix = fns.str();
		model.suffix_serial = serial.str();
	}

	void setup_sys_size(Model& model) override
	{
		const int N = model.ini.GetInteger("lind_ham", "N", 0);
		model.sys_size = N;
	}

	void setup_period(Model& model) override
	{
	}

	void setup_hamiltonian(Model& model) override
	{
		model.hamiltonian = get_H_mtx(model, model.sys_size);
	}

	void setup_hamiltonian_drv(Model& model) override
	{
		model.log_message("hamiltonian_drv is absent in this model");
	}

	void setup_dissipators(Model& model) override
	{
		int M = model.sys_size * model.sys_size - 1;

		ds_mtx G = get_G_mtx(model, model.sys_size * model.sys_size - 1, model.sys_size * model.sys_size - 1);

		model.log_time_duration();
		model.log_message("G eigen...");

		Eigen::ComplexEigenSolver<ds_mtx> es;
		es.compute(G, true);

		auto evals_tmp = es.eigenvalues();
		std::vector<std::complex<double>> G_evals(evals_tmp.data(), evals_tmp.data() + evals_tmp.rows() * evals_tmp.cols());
		auto evecs = es.eigenvectors();

		for (int k1 = 0; k1 < M; k1++)
		{
			sp_mtx diss;
			for (int k2 = 0; k2 < M; k2++)
			{
				diss += evecs(k2, k1) * model.f_basis[k2 + 1];
			}
			diss *= std::sqrt(evals_tmp[k1]);
			
			model.dissipators.push_back(diss);
		}
	}

	void setup_lindbladian(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		const auto alpha = model.ini.GetReal("lind_ham", "alpha", 0.0);

		const std::complex<double> i1(0.0, 1.0);
		const sp_mtx eye = get_sp_eye(model.sys_size);
		const sp_mtx hamiltonian_transposed(model.hamiltonian.transpose());

		model.lindbladian = -i1 * (Eigen::kroneckerProduct(eye, model.hamiltonian) - Eigen::kroneckerProduct(hamiltonian_transposed, eye)) * alpha / std::sqrt(static_cast<double>(model.sys_size));

		for (const auto& diss : model.dissipators)
		{
			sp_mtx diss_tmp_1((diss.adjoint()).transpose());
			sp_mtx diss_tmp_2(diss.adjoint() * diss);
			sp_mtx diss_tmp_3(diss_tmp_2.transpose());

			model.lindbladian += 0.5 * (2.0 *
				Eigen::kroneckerProduct(eye, diss) *
				Eigen::kroneckerProduct(diss_tmp_1, eye) -
				Eigen::kroneckerProduct(diss_tmp_3, eye) -
				Eigen::kroneckerProduct(eye, diss_tmp_2));
		}
	}

	void setup_lindbladians_drv(Model& model) override
	{
		model.log_message("lindbladians_drv is absent in this model");
	}

	void release_observables(Model& model) override
	{
	}

	static ds_mtx get_G_mtx(Model& model, size_t dim_main, size_t dim_aux)
	{
		const std::string run_type = model.ini.Get("global", "run_type", "unknown");

		int seed;

		if (run_type == "serial")
		{
			seed = std::round(model.serial_state);
		}
		else
		{
			seed = model.ini.GetInteger("lind_ham", "seed", 0);
		}

		const int num_seeds = model.ini.GetInteger("lind_ham", "num_seeds", 0);

		VSLStreamStatePtr stream;
		vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
		vslLeapfrogStream(stream, seed, num_seeds);
		double* disorder_real = new double[dim_main * dim_aux];
		double* disorder_imag = new double[dim_main * dim_aux];
		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, dim_main * dim_aux, disorder_real, 0.0, 1.0);
		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, dim_main * dim_aux, disorder_imag, 0.0, 1.0);
		model.log_message("G random complete");
		model.log_time_duration();

		ds_mtx X(dim_main, dim_aux);

		for (auto st_id_1 = 0; st_id_1 < dim_main; st_id_1++)
		{
			for (auto st_id_2 = 0; st_id_2 < dim_aux; st_id_2++)
			{
				auto index = st_id_1 * dim_aux + st_id_2;
				X(st_id_1, st_id_2) = std::complex<double>(0.5 * disorder_real[index], 0.5 * disorder_imag[index]);
			}
		}

		delete[] disorder_real;
		delete[] disorder_imag;

		ds_mtx G = X * X.adjoint();
		model.log_message(fmt::format("G_size = {:d} rows,  {:d} cols ", G.rows(), G.cols()));
		const auto trace = G.trace();
		G = double(model.sys_size) * G / trace.real();

		model.log_message("G generation complete");
		model.log_time_duration();

		return G;
	}

	static ds_mtx get_H_mtx(Model& model, size_t dim)
	{
		const std::string run_type = model.ini.Get("global", "run_type", "unknown");

		int seed;

		if (run_type == "serial")
		{
			seed = std::round(model.serial_state);
		}
		else
		{
			seed = model.ini.GetInteger("lind_ham", "seed", 0);
		}

		const int num_seeds = model.ini.GetInteger("lind_ham", "num_seeds", 0);

		VSLStreamStatePtr stream;
		vslNewStream(&stream, VSL_BRNG_MCG31, 13371337);
		vslLeapfrogStream(stream, seed, num_seeds);
		double* disorder = new double[dim * dim];
		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, dim * dim, disorder, 0.0, 1.0);
		ds_mtx X(dim, dim);
		for (auto st_id_1 = 0; st_id_1 < dim; st_id_1++)
		{
			for (auto st_id_2 = 0; st_id_2 < dim; st_id_2++)
			{
				auto index = st_id_1 * dim + st_id_2;
				X(st_id_1, st_id_2) = std::complex<double>(disorder[index], disorder[index]);
			}
		}
		delete[] disorder;

		ds_mtx y = 0.5 * (X + X.adjoint());
		ds_mtx y2 = y * y;
		const auto y2_tr = y2.trace();
		y = y / std::sqrt(y2_tr.real());

		return y;
	}

	void setup_serial_data(
		Model& model,
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex) override
	{
		features_complex.insert({ "rho_evals", {} });
		features_complex.insert({ "lindbladian_evals", {} });
	}

	void fill_serial_data(
		Model& model,
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex) override
	{
		features_complex["rho_evals"].insert(features_complex["rho_evals"].end(), model.rho_evals.begin(), model.rho_evals.end());
		features_complex["lindbladian_evals"].insert(features_complex["lindbladian_evals"].end(), model.lindbladian_evals.begin(), model.lindbladian_evals.end());
	}
};
