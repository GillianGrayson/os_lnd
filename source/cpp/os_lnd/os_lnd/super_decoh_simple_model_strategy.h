#pragma once
#include "init.h"
#include "model.h"
#include <gcem.hpp>
#include <mkl.h>
#include <iomanip>
#include "save.h"
#include "model_strategy.h"
#include <unsupported/Eigen/KroneckerProduct>
#include "routines.h"
#include <numeric>

struct SuperDecohSimpleModelStrategy : ModelStrategy
{
	void setup_aux_data(Model& model) override
	{
	}

	void setup_suffix(Model& model) override
	{
		const int name_precision = model.ini.GetInteger("global", "name_precision", 0);

		const int N = model.ini.GetInteger("super_decoh_simple", "N", 0);
		const auto p = model.ini.GetReal("super_decoh_simple", "p", 0.0);

		const int seed = model.ini.GetInteger("super_decoh_simple", "seed", 0);

		const int reshuffle_type = model.ini.GetInteger("super_decoh_simple", "reshuffle_type", 0);

		std::stringstream fns;
		fns << "_reshuffle(" << reshuffle_type << ")";
		fns << "_N(" << N << ")";
		fns << "_p(" << std::setprecision(name_precision) << std::fixed << p << ")";
		fns << "_seed(" << seed << ")";
		fns << ".txt";

		model.suffix = fns.str();
	}

	void setup_sys_size(Model& model) override
	{
		const int N = model.ini.GetInteger("super_decoh_simple", "N", 0);
		model.sys_size = N;
	}

	void setup_period(Model& model) override
	{
	}

	void setup_hamiltonian(Model& model) override
	{
		model.log_message("hamiltonian is absent in this model");
	}

	void setup_hamiltonian_drv(Model& model) override
	{
		model.log_message("hamiltonian_drv is absent in this model");
	}

	void setup_dissipators(Model& model) override
	{
		model.log_message("dissipators are absent in this model");
	}

	void setup_lindbladian(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		const int reshuffle_type = model.ini.GetInteger("super_decoh_simple", "reshuffle_type", 0);

		ds_mtx rho_tmp = get_rho_mtx(model);

		if (debug_dump || save_G)
		{
			auto fn = "G_mtx" + model.suffix;
			save_dense_mtx(G, fn, save_precision);
		}

		const sp_mtx eye = get_sp_eye(model.sys_size);

		model.lindbladian_dense = ds_mtx::Zero(model.sys_size * model.sys_size, model.sys_size * model.sys_size);

		auto M = model.sys_size * model.sys_size - 1;
		sp_mtx tmp;

		std::vector<sp_mtx> lefts(M);
		std::vector<sp_mtx> rights(M);
		for (auto k1 = 0; k1 < M; k1++)
		{
			lefts[k1] = Eigen::kroneckerProduct(model.f_basis[k1 + 1], eye);
			rights[k1] = Eigen::kroneckerProduct(eye, model.f_basis[k1 + 1].transpose());
		}
		model.log_message("Lefts and rights created");
		model.log_time_duration();

		for (auto k1 = 0; k1 < M; k1++)
		{
			for (auto k2 = 0; k2 < M; k2++)
			{
				model.lindbladian_dense += G(k1, k2) * lefts[k1] * rights[k2];
			}
		}

		ds_mtx reshuffle;
		if (reshuffle_type == 1)
		{
			reshuffle = get_reshuffle_ds_mtx_1(model.lindbladian_dense, model.sys_size * model.sys_size, model.sys_size);
		}
		else if (reshuffle_type == 0)
		{
			reshuffle = get_reshuffle_ds_mtx_0(model.lindbladian_dense, model.sys_size * model.sys_size, model.sys_size);
		}
		else
		{
			model.throw_error("Unsupported reshuffle_type");
		}
		model.lindbladian_dense = reshuffle;

		decoherence(model, model.lindbladian_dense);

		ds_mtx A = get_addition_ds_mtx(model.lindbladian_dense, model.sys_size * model.sys_size, model.sys_size);

		if (debug_dump || save_A)
		{
			auto fn = "A_mtx" + model.suffix;
			save_dense_mtx(A, fn, save_precision);
		}

		if (reshuffle_type == 1)
		{
			reshuffle = get_reshuffle_ds_mtx_1(model.lindbladian_dense, model.sys_size * model.sys_size, model.sys_size);
		}
		else if (reshuffle_type == 0)
		{
			reshuffle = get_reshuffle_ds_mtx_0(model.lindbladian_dense, model.sys_size * model.sys_size, model.sys_size);
		}
		else
		{
			model.throw_error("Unsupported reshuffle_type");
		}
		model.lindbladian_dense = reshuffle;

		model.lindbladian_dense -= 0.5 * (Eigen::kroneckerProduct(A, eye) + Eigen::kroneckerProduct(eye, A.transpose()));
	}

	void setup_lindbladian_drv(Model& model) override
	{
		model.log_message("lindbladian_drv is absent in this model");
	}

	void release_observables(Model& model) override
	{
	}

	static ds_mtx get_rho_mtx(Model& model)
	{
		const int seed = model.ini.GetInteger("super_decoh_simple", "seed", 0);
		const int num_seeds = model.ini.GetInteger("super_decoh_simple", "num_seeds", 0);

		int M = model.sys_size * model.sys_size - 1;

		VSLStreamStatePtr stream;
		vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
		vslLeapfrogStream(stream, seed, num_seeds);

		double* disorder_real = new double[M * M];
		double* disorder_imag = new double[M * M];

		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, M * M, disorder_real, 0.0, 1.0);
		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, M * M, disorder_imag, 0.0, 1.0);

		model.log_message("rho random complete");
		model.log_time_duration();

		ds_mtx X(M, M);

		for (auto st_id_1 = 0; st_id_1 < M; st_id_1++)
		{
			for (auto st_id_2 = 0; st_id_2 < M; st_id_2++)
			{
				auto index = st_id_1 * M + st_id_2;
				X(st_id_1, st_id_2) = std::complex<double>(0.5 * disorder_real[index], 0.5 * disorder_imag[index]);
			}
		}

		delete[] disorder_real;
		delete[] disorder_imag;

		ds_mtx rho = X * X.adjoint();

		const auto trace = rho.trace();
		rho = double(model.sys_size) * rho / trace.real();
		model.lindbladian_evals_mult = std::complex<double>(1.0, 0.0);

		model.log_message("rho generation complete");
		model.log_time_duration();

		return rho;
	}

	static void decoherence(Model& model, ds_mtx& mtx)
	{
		const auto p = model.ini.GetReal("super_decoh", "p", 0.0);

		Eigen::VectorXcd diag_vec(mtx.diagonal());

		mtx = p * mtx;

		mtx.diagonal() = diag_vec;
	}
};
