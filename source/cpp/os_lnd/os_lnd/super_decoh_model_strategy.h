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

		const int reshuffle_type = model.ini.GetInteger("super_decoh", "reshuffle_type", 0);
		const int G_type = model.ini.GetInteger("super_decoh", "G_type", 0);

		const int aux_dim = model.ini.GetInteger("super_decoh", "aux_dim", 0);

		std::stringstream fns;
		fns << "_reshuffle(" << reshuffle_type << ")";
		fns << "_G(" << G_type << ")";
		fns << "_N(" << N << ")";
		fns << "_ad(" << aux_dim << ")";
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
		const int G_type = model.ini.GetInteger("super_decoh", "G_type", 0);
		const auto save_G = model.ini.GetBoolean("super_decoh", "save_G", false);
		const auto save_A = model.ini.GetBoolean("super_decoh", "save_A", false);
		const int reshuffle_type = model.ini.GetInteger("super_decoh", "reshuffle_type", 0);
		std::string method = model.ini.Get("super_decoh", "method", "origin");

		if (method == "origin")
		{
			model.init_f_basis();
			
			ds_mtx G;
			if (G_type == 1)
			{
				const int aux_dim = model.ini.GetInteger("super_decoh", "aux_dim", 0);
				G = get_G_mtx(model, model.sys_size * model.sys_size - 1, aux_dim);
			}
			else if (G_type == 0)
			{
				G = get_G_mtx(model, model.sys_size * model.sys_size - 1, model.sys_size * model.sys_size - 1);
			}
			else
			{
				model.throw_error("Unsupported G_type");
			}

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
		else if (method == "simple")
		{
			ds_mtx G = get_G_mtx(model, model.sys_size * model.sys_size, model.sys_size * model.sys_size);

			if (debug_dump || save_G)
			{
				auto fn = "G_mtx" + model.suffix;
				save_dense_mtx(G, fn, save_precision);
			}

			decoherence(model, G);

			ds_mtx reshuffle;
			if (reshuffle_type == 1)
			{
				reshuffle = get_reshuffle_ds_mtx_1(G, model.sys_size * model.sys_size, model.sys_size);
			}
			else if (reshuffle_type == 0)
			{
				reshuffle = get_reshuffle_ds_mtx_0(G, model.sys_size * model.sys_size, model.sys_size);
			}
			else
			{
				model.throw_error("Unsupported reshuffle_type");
			}

			ds_mtx A = ds_mtx::Zero(model.sys_size, model.sys_size);
			for (int st_1 = 0; st_1 < model.sys_size; st_1++)
			{
				for (int st_2 = 0; st_2 < model.sys_size; st_2++)
				{
					for (int st_3 = 0; st_3 < model.sys_size; st_3++)
					{
						int index_1 = st_3 + model.sys_size * st_1;
						int index_2 = st_3 + model.sys_size * st_2;

						A(st_1, st_2) += G(index_1, index_2);
					}
				}
			}
			
			ds_mtx eye = ds_mtx::Identity(model.sys_size, model.sys_size);
			
			model.lindbladian_dense = reshuffle;
			model.lindbladian_dense -= 0.5 * (Eigen::kroneckerProduct(A, eye) + Eigen::kroneckerProduct(eye, A.transpose()));
		}
		else
		{
			model.throw_error("Unsupported method");
		}
	}

	void setup_lindbladian_drv(Model& model) override
	{
		model.log_message("lindbladian_drv is absent in this model");
	}

	void release_observables(Model& model) override
	{
	}

	static ds_mtx get_G_mtx(Model& model, size_t dim_main, size_t dim_aux)
	{
		const int seed = model.ini.GetInteger("super_decoh", "seed", 0);
		const int num_seeds = model.ini.GetInteger("super_decoh", "num_seeds", 0);
		const auto evals_G = model.ini.GetBoolean("super_decoh", "evals_G", false);
		
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

		if (evals_G)
		{
			if (dim_main == dim_aux)
			{
				calc_evals_G(model, X);
			}
			else
			{
				model.log_message("Can't calc_evals_G");
			}
		}

		ds_mtx G = X * X.adjoint();

		model.log_message(fmt::format("G_size = {:d} rows,  {:d} cols ", G.rows(), G.cols()));

		const auto trace = G.trace();
		G = double(model.sys_size) * G / trace.real();

		model.log_message("G generation complete");
		model.log_time_duration();

		return G;
	}

	static void decoherence(Model& model, ds_mtx& mtx)
	{
		const auto p = model.ini.GetReal("super_decoh", "p", 0.0);

		Eigen::VectorXcd diag_vec(mtx.diagonal());

		mtx = p * mtx;

		mtx.diagonal() = diag_vec;
	}

	static void calc_evals_G(Model& model, ds_mtx& mtx)
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		
		model.log_time_duration();
		model.log_message("G eigen...");

		Eigen::ComplexEigenSolver<ds_mtx> es;
		es.compute(mtx, true);

		auto evals_tmp = es.eigenvalues();
		std::vector<std::complex<double>> evals(evals_tmp.data(), evals_tmp.data() + evals_tmp.rows() * evals_tmp.cols());
		auto fn = "G_evals" + model.suffix;
		save_vector(evals, fn, save_precision);

		model.log_time_duration();
		model.log_message("G eigen done");
	}
};
