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

struct SuperDecohModelStrategy : ModelStrategy
{
	void setup_aux_data(Model& model) override
	{
	}

	void setup_suffix(Model& model) override
	{
		const int name_precision = model.ini.GetInteger("global", "name_precision", 0);

		const int N = model.ini.GetInteger("super_decoh", "N", 0);
		const auto p = model.ini.GetReal("super_decoh", "p", 0.0);

		const int seed = model.ini.GetInteger("super_decoh", "seed", 0);

		std::stringstream fns;
		fns << "_N(" << N << ")";
		fns << "_p(" << std::setprecision(name_precision) << std::fixed << p << ")";
		fns << "_seed(" << seed << ")";
		fns << ".txt";

		model.suffix = fns.str();
	}

	void setup_sys_size(Model& model) override
	{
		const int N = model.ini.GetInteger("super_decoh", "N", 0);
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
		const auto save_G = model.ini.GetBoolean("super_decoh", "save_G", false);
		
		model.init_f_basis();

		ds_mtx G = get_G_mtx(model);

		if (debug_dump || save_G)
		{
			auto fn = "G_mtx" + model.suffix;
			save_dense_mtx(G, fn, save_precision);
		}
		
		const sp_mtx eye = get_sp_eye(model.sys_size);

		model.lindbladian_dense = ds_mtx::Zero(model.sys_size * model.sys_size, model.sys_size * model.sys_size);

		auto M = model.sys_size * model.sys_size - 1;
		sp_mtx tmp;
		for (auto k1 = 0; k1 < M; k1++)
		{
			for (auto k2 = 0; k2 < M; k2++)
			{
				tmp = G(k1, k2) * Eigen::kroneckerProduct(model.f_basis[k1 + 1], eye) * Eigen::kroneckerProduct(eye, model.f_basis[k2 + 1].transpose());
				model.lindbladian_dense += tmp;
			}
		}

		ds_mtx reshuffle = get_reshuffle_ds_mtx(model.lindbladian_dense, model.sys_size * model.sys_size, model.sys_size);
		model.lindbladian_dense = reshuffle;
	}

	void setup_lindbladian_drv(Model& model) override
	{
		model.log_message("lindbladian_drv is absent in this model");
	}

	void release_observables(Model& model) override
	{
	}

	static ds_mtx get_G_mtx(Model& model)
	{
		const int seed = model.ini.GetInteger("super_decoh", "seed", 0);
		const int num_seeds = model.ini.GetInteger("super_decoh", "num_seeds", 0);

		int M = model.sys_size * model.sys_size - 1;
		
		VSLStreamStatePtr stream;
		vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
		vslLeapfrogStream(stream, seed, num_seeds);
		
		double* disorder_real = new double[M * M];
		double* disorder_imag = new double[M * M];

		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, M * M, disorder_real, 0.0, 1.0);
		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, M * M, disorder_imag, 0.0, 1.0);

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

		ds_mtx G = X * X.adjoint();

		const auto trace = G.trace();

		G = double(model.sys_size) * G / trace.real();

		return G;
	}
};
