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
#include <Eigen/SVD>
#include "routines.h"
#include <numeric>


struct XXZModelStrategy : ModelStrategy
{
	std::vector<double> energies;
	std::vector<sp_mtx> sigma_x_mtxs;
	std::vector<sp_mtx> sigma_y_mtxs;
	std::vector<sp_mtx> sigma_z_mtxs;
	std::vector<sp_mtx> sigma_m_mtxs;
	std::vector<sp_mtx> sigma_p_mtxs;
	sp_mtx jznd_mtx;
	sp_mtx jvak_mtx;

	void setup_aux_data(Model& model) override
	{
		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		
		const int quantity_index = model.ini.GetInteger("xxz", "quantity_index", 0);
		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);

		const std::complex<double> i1(0.0, 1.0);

		sigma_x_mtxs = get_kronecker_mtxs(num_spins, "sigma_x");
		sigma_y_mtxs = get_kronecker_mtxs(num_spins, "sigma_y");
		sigma_z_mtxs = get_kronecker_mtxs(num_spins, "sigma_z");
		sigma_m_mtxs = get_kronecker_mtxs(num_spins, "sigma_m");
		sigma_p_mtxs = get_kronecker_mtxs(num_spins, "sigma_p");

		jznd_mtx = (sigma_x_mtxs[quantity_index] * sigma_y_mtxs[quantity_index + 1] - sigma_y_mtxs[quantity_index] * sigma_x_mtxs[quantity_index + 1]);

		jvak_mtx = sp_mtx(model.sys_size, model.sys_size);
		for (auto spin_id = 0; spin_id < num_spins - 1; spin_id++)
		{
			jvak_mtx += (sigma_p_mtxs[spin_id] * sigma_m_mtxs[spin_id + 1] - sigma_m_mtxs[spin_id] * sigma_p_mtxs[spin_id + 1]);
		}
		jvak_mtx = 2.0 * i1 / double(num_spins - 1) * jvak_mtx;

		if (debug_dump)
		{
			auto fn = "jznd_mtx" + model.suffix;
			save_sp_mtx(jznd_mtx, fn, save_precision);

			fn = "jvak_mtx" + model.suffix;
			save_sp_mtx(jvak_mtx, fn, save_precision);
		}
	}

	void setup_suffix(Model& model) override
	{
		const int name_precision = model.ini.GetInteger("global", "name_precision", 0);

		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);

		const int seed = model.ini.GetInteger("xxz", "seed", 0);

		const auto W = model.ini.GetReal("xxz", "W", 0.0);
		const auto mu = model.ini.GetReal("xxz", "mu", 0.0);
		const int drv_type = model.ini.GetInteger("xxz", "drv_type", 0);
		const auto T1 = model.ini.GetReal("xxz", "T1", 0.0);
		const auto T2 = model.ini.GetReal("xxz", "T2", 0.0);

		const int quantity_index = model.ini.GetInteger("xxz", "quantity_index", 0);

		std::stringstream fns;
		fns << "_ns(" << num_spins << ")";
		fns << "_seed(" << seed << ")";
		fns << "_prm(" << std::setprecision(name_precision) << std::fixed << W << "_" << mu << "_" << drv_type << "_" << T1 << "_" << T2 << ")";
		fns << "_j(" << quantity_index << ")";
		fns << ".txt";

		model.suffix = fns.str();
	}

	void setup_sys_size(Model& model) override
	{
		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);

		model.sys_size = std::pow(2, num_spins);
	}

	void setup_period(Model& model) override
	{
		auto T1 = model.ini.GetReal("xxz", "T1", 0.0);
		auto T2 = model.ini.GetReal("xxz", "T2", 0.0);

		model.period = std::max(T1, T2);
	}

	void setup_hamiltonian(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		
		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);
		
		const int seed = model.ini.GetInteger("xxz", "seed", 0);
		const int num_seeds = model.ini.GetInteger("xxz", "num_seeds", 0);
		
		const auto W = model.ini.GetReal("xxz", "W", 0.0);

		model.hamiltonian = sp_mtx(model.sys_size, model.sys_size);

		energies = std::vector<double>(num_spins, 0.0);

		VSLStreamStatePtr stream;
		vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
		vslLeapfrogStream(stream, seed, num_seeds);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, num_spins, energies.data(), -1.0, 1.0);

		auto fn = "energies" + model.suffix;
		save_vector(energies, fn, save_precision);

		for (auto spin_id = 0; spin_id < num_spins; spin_id++)
		{
			if (spin_id < num_spins - 1)
			{
				model.hamiltonian += (sigma_x_mtxs[spin_id] * sigma_x_mtxs[spin_id + 1] + sigma_y_mtxs[spin_id] * sigma_y_mtxs[spin_id + 1] + sigma_z_mtxs[spin_id] * sigma_z_mtxs[spin_id + 1]);
			}
			model.hamiltonian += W * energies[spin_id] * sigma_z_mtxs[spin_id];
		}
	}

	void setup_hamiltonian_drv(Model& model) override
	{
		model.log_message("hamiltonian_drv is absent in this model");
	}

	void setup_dissipators(Model& model) override
	{
		const int num_spins = model.ini.GetInteger("mbl", "num_spins", 0);
		
		int num_diss = 4;
		model.dissipators.push_back(sigma_p_mtxs[0]);
		model.dissipators.push_back(sigma_m_mtxs[0]);
		model.dissipators.push_back(sigma_p_mtxs[num_spins - 1]);
		model.dissipators.push_back(sigma_m_mtxs[num_spins - 1]);
	}

	void setup_lindbladian(Model& model) override
	{
		const std::complex<double> i1(0.0, 1.0);
		const sp_mtx eye = get_sp_eye(model.sys_size);
		const sp_mtx hamiltonian_transposed(model.hamiltonian.transpose());

		model.lindbladian = -i1 * (Eigen::kroneckerProduct(eye, model.hamiltonian) - Eigen::kroneckerProduct(hamiltonian_transposed, eye));

		for (auto diss_id = 0; diss_id != model.dissipators.size(); diss_id++) 
		{
			sp_mtx diss = model.dissipators[diss_id];

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
		const sp_mtx eye = get_sp_eye(model.sys_size);

		for (auto diss_id = 0; diss_id != model.dissipators.size(); diss_id++)
		{
			sp_mtx diss = model.dissipators[diss_id];

			sp_mtx diss_tmp_1((diss.adjoint()).transpose());
			sp_mtx diss_tmp_2(diss.adjoint() * diss);
			sp_mtx diss_tmp_3(diss_tmp_2.transpose());

			sp_mtx tmp = 0.5 * (2.0 *
				Eigen::kroneckerProduct(eye, diss) *
				Eigen::kroneckerProduct(diss_tmp_1, eye) -
				Eigen::kroneckerProduct(diss_tmp_3, eye) -
				Eigen::kroneckerProduct(eye, diss_tmp_2));

			model.lindbladians_drv.push_back(tmp);
		}
	}

	double get_quantity_znd(Model& model)
	{
		Eigen::MatrixXcd op = jznd_mtx * model.rho;
		double quantity = op.trace().real();
		return quantity;
	}

	double get_quantity_vak(Model& model)
	{
		Eigen::MatrixXcd op = jvak_mtx * model.rho;
		double quantity = op.trace().real();
		return quantity;
	}

	void release_observables(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		std::string fn;

		const double jznd = get_quantity_znd(model);
		fn = "znd" + model.suffix;
		save_value(jznd, fn, save_precision);

		const double jvak = get_quantity_vak(model);
		fn = "vak" + model.suffix;
		save_value(jvak, fn, save_precision);
	}
};
