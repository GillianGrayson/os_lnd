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

struct LFKModelStrategy : ModelStrategy
{
	void setup_aux_data(Model& model) override
	{
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
			seed = model.ini.GetInteger("lfk", "seed", 0);
		}

		const int name_precision = model.ini.GetInteger("global", "name_precision", 0);

		const int N = model.ini.GetInteger("lfk", "N", 0);

		std::stringstream fns;
		fns << "_N(" << N << ")";

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
		const int N = model.ini.GetInteger("lfk", "N", 0);
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

		const std::string run_type = model.ini.Get("global", "run_type", "unknown");
		int seed;
		if (run_type == "serial")
		{
			seed = std::round(model.serial_state);
		}
		else
		{
			seed = model.ini.GetInteger("lfk", "seed", 0);
		}

		const int num_seeds = model.ini.GetInteger("lfk", "num_seeds", 0);

		VSLStreamStatePtr stream;
		vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
		vslLeapfrogStream(stream, seed, num_seeds);

		double* MR = new double[model.sys_size * model.sys_size];
		double* MI = new double[model.sys_size * model.sys_size];

		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, model.sys_size * model.sys_size, MR, 0.0, 1.0);
		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, model.sys_size * model.sys_size, MI, 0.0, 1.0);

		ds_mtx MS(model.sys_size, model.sys_size);
		for (auto st_id_1 = 0; st_id_1 < model.sys_size; st_id_1++)
		{
			for (auto st_id_2 = 0; st_id_2 < model.sys_size; st_id_2++)
			{
				auto index = st_id_1 * model.sys_size + st_id_2;
				MS(st_id_1, st_id_2) = std::complex<double>(MR[index], MI[index]);
			}
		}
		delete[] MR;
		delete[] MI;

		if (debug_dump)
		{
			auto fn = "MS_mtx" + model.suffix;
			save_dense_mtx(MS, fn, save_precision);
		}
		
		ds_mtx MZ(model.sys_size, model.sys_size);
		for (auto st_id_1 = 0; st_id_1 < model.sys_size; st_id_1++)
		{
			for (auto st_id_2 = 0; st_id_2 < model.sys_size; st_id_2++)
			{
				auto index = st_id_1 * model.sys_size + st_id_2;
				MZ(st_id_2, st_id_1) = MS(st_id_2, st_id_1) * std::conj(MS(st_id_2, st_id_1));
			}
		}
		auto trace_MZ = MZ.trace();
		MZ = std::sqrt(static_cast<double>(model.sys_size)) * MZ / trace_MZ;

		double* disorder_real = new double[model.sys_size * model.sys_size * model.sys_size * model.sys_size];
		double* disorder_imag = new double[model.sys_size * model.sys_size * model.sys_size * model.sys_size];

		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, model.sys_size * model.sys_size * model.sys_size * model.sys_size, disorder_real, 0.0, 1.0);
		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, model.sys_size * model.sys_size * model.sys_size * model.sys_size, disorder_imag, 0.0, 1.0);

		ds_mtx X(model.sys_size * model.sys_size, model.sys_size * model.sys_size);

		for (auto st_id_1 = 0; st_id_1 < model.sys_size * model.sys_size; st_id_1++)
		{
			for (auto st_id_2 = 0; st_id_2 < model.sys_size * model.sys_size; st_id_2++)
			{
				auto index = st_id_1 * model.sys_size * model.sys_size + st_id_2;
				X(st_id_1, st_id_2) = std::complex<double>(disorder_real[index], disorder_imag[index]);
			}
		}
		delete[] disorder_real;
		delete[] disorder_imag;

		if (debug_dump)
		{
			auto fn = "X_mtx" + model.suffix;
			save_dense_mtx(X, fn, save_precision);
		}

		ds_mtx rho = X * X.adjoint();

		for (auto j = 0; j < model.sys_size; j++)
		{
			for (auto k = 0; k < model.sys_size; k++)
			{
				auto s = model.sys_size * j + k;
				auto ttt = rho(s, s) / MZ(j, k);

				for (auto sp = 0; sp < model.sys_size * model.sys_size; sp++)
				{
					X(sp, s) = X(sp, s) / std::sqrt(ttt);
				}
			}
		}

		rho = X * X.adjoint();
		auto trace_rho = rho.trace();
		rho = std::sqrt(static_cast<double>(model.sys_size * model.sys_size)) * rho / trace_rho;

		model.lindbladian_dense = get_reshuffle_ds_mtx_0(rho, model.sys_size * model.sys_size, model.sys_size);

		ds_mtx AA = ds_mtx::Zero(model.sys_size, model.sys_size);
		for (auto s1 = 0; s1 < model.sys_size; s1++)
		{
			for (auto s2 = 0; s2 < model.sys_size; s2++)
			{
				for (auto s3 = 0; s3 < model.sys_size; s3++)
				{
					auto w1 = s3 + model.sys_size * s1;
					auto w2 = s3 + model.sys_size * s2;
					AA(s1, s2) = AA(s1, s2) + rho(w1, w2);
				}
			}
		}
		const sp_mtx eye = get_sp_eye(model.sys_size);
		model.lindbladian_dense -= 0.5 * (Eigen::kroneckerProduct(AA, eye) + Eigen::kroneckerProduct(eye, AA.transpose()));
	}

	void setup_lindbladians_drv(Model& model) override
	{
		model.log_message("lindbladians_drv is absent in this model");
	}

	void release_observables(Model& model) override
	{
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
