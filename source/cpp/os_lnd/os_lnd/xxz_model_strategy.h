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

	std::vector<sp_mtx> j_mtxs;

	void setup_aux_data(Model& model) override
	{
		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		
		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);

		const std::complex<double> i1(0.0, 1.0);

		sigma_x_mtxs = get_kronecker_mtxs(num_spins, "sigma_x");
		sigma_y_mtxs = get_kronecker_mtxs(num_spins, "sigma_y");
		sigma_z_mtxs = get_kronecker_mtxs(num_spins, "sigma_z");
		sigma_m_mtxs = get_kronecker_mtxs(num_spins, "sigma_m");
		sigma_p_mtxs = get_kronecker_mtxs(num_spins, "sigma_p");

		for (auto spin_id = 0; spin_id < num_spins - 1; spin_id++)
		{
			sp_mtx j_mtx = 0.25 * (sigma_x_mtxs[spin_id] * sigma_y_mtxs[spin_id + 1] - sigma_y_mtxs[spin_id] * sigma_x_mtxs[spin_id + 1]);
			j_mtxs.push_back(j_mtx);
			if (debug_dump)
			{
				auto fn = "znd_mtx_" + std::to_string(spin_id) + model.suffix;
				save_sp_mtx(j_mtx, fn, save_precision);
			}
		}


		const std::string run_type = model.ini.Get("global", "run_type", "unknown");
		int seed;
		if (run_type == "serial")
		{
			seed = std::round(model.serial_state);
		}
		else
		{
			seed = model.ini.GetInteger("xxz", "seed", 0);
		}
		const int num_seeds = model.ini.GetInteger("xxz", "num_seeds", 0);

		energies = std::vector<double>(num_spins, 0.0);

		VSLStreamStatePtr stream;
		vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
		vslLeapfrogStream(stream, seed, num_seeds);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, num_spins, energies.data(), -1.0, 1.0);

		if (debug_dump)
		{
			auto fn = "energies" + model.suffix;
			save_vector(energies, fn, save_precision);
		}
	}

	void setup_suffix(Model& model) override
	{
		const std::string run_type = model.ini.Get("global", "run_type", "unknown");
		const int name_precision = model.ini.GetInteger("global", "name_precision", 0);
		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);

		int seed;
		if (run_type == "serial")
		{
			seed = std::round(model.serial_state);
		}
		else
		{
			seed = model.ini.GetInteger("xxz", "seed", 0);
		}

		const auto W = model.ini.GetReal("xxz", "W", 0.0);
		const auto Delta = model.ini.GetReal("xxz", "Delta", 0.0);
		const auto mu = model.ini.GetReal("xxz", "mu", 0.0);
		const int drv_type = model.ini.GetInteger("xxz", "drv_type", 0);
		const auto ampl = model.ini.GetReal("xxz", "ampl", 0.0);
		const auto freq = model.ini.GetReal("xxz", "freq", 0.0);
		const auto phase = model.ini.GetReal("xxz", "phase", 0.0);

		std::stringstream fns;
		fns << "_ns(" << num_spins << ")";
		fns << "_prm(" << std::setprecision(name_precision) << std::fixed << Delta << "_" << W << "_" << mu << "_" << drv_type << "_" << ampl << "_" << freq << "_" << phase << ")";

		std::stringstream serial;
		serial << fns.rdbuf();
		fns << "_seed(" << seed << ")";

		fns << ".txt";
		serial << ".txt";

		model.suffix = fns.str();
		model.suffix_serial = serial.str();
	}

	void setup_sys_size(Model& model) override
	{
		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);
		model.sys_size = std::pow(2, num_spins);
	}

	void setup_period(Model& model) override
	{
		const auto freq = model.ini.GetReal("xxz", "freq", 0.0);
		double pi = std::atan(1.0) * 4.0;
		model.period = 2.0 * pi / freq;
	}

	void setup_hamiltonian(Model& model) override
	{
		const std::string run_type = model.ini.Get("global", "run_type", "unknown");
		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		
		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);
		
		int seed;
		if (run_type == "serial")
		{
			seed = std::round(model.serial_state);
		}
		else
		{
			seed = model.ini.GetInteger("xxz", "seed", 0);
		}
		
		const int num_seeds = model.ini.GetInteger("xxz", "num_seeds", 0);

		const auto Delta = model.ini.GetReal("xxz", "Delta", 0.0);
		const auto W = model.ini.GetReal("xxz", "W", 0.0);

		model.hamiltonian = sp_mtx(model.sys_size, model.sys_size);

		energies = std::vector<double>(num_spins, 0.0);

		VSLStreamStatePtr stream;
		vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
		vslLeapfrogStream(stream, seed, num_seeds);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, num_spins, energies.data(), -1.0, 1.0);

		if (debug_dump)
		{
			auto fn = "energies" + model.suffix;
			save_vector(energies, fn, save_precision);
		}

		for (auto spin_id = 0; spin_id < num_spins - 1; spin_id++)
		{
			model.hamiltonian += (0.25 * sigma_x_mtxs[spin_id] * sigma_x_mtxs[spin_id + 1] + 0.25 * sigma_y_mtxs[spin_id] * sigma_y_mtxs[spin_id + 1] + 0.25 * Delta * sigma_z_mtxs[spin_id] * sigma_z_mtxs[spin_id + 1] + 0.25 * W * energies[spin_id] * sigma_z_mtxs[spin_id] + 0.25 * W * energies[spin_id + 1] * sigma_z_mtxs[spin_id + 1]);
		}
	}

	void setup_hamiltonian_drv(Model& model) override
	{
		model.log_message("hamiltonian_drv is absent in this model");
	}

	void setup_dissipators(Model& model) override
	{
		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);
		
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

			model.lindbladian += 0.25 * (2.0 *
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

			sp_mtx tmp = 0.25 * (2.0 *
				Eigen::kroneckerProduct(eye, diss) *
				Eigen::kroneckerProduct(diss_tmp_1, eye) -
				Eigen::kroneckerProduct(diss_tmp_3, eye) -
				Eigen::kroneckerProduct(eye, diss_tmp_2));

			model.lindbladians_drv.push_back(tmp);
		}
	}

	std::vector<double> get_quantities(Model& model)
	{
		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);

		std::vector<double> quantities;
		for (auto spin_id = 0; spin_id < num_spins - 1; spin_id++)
		{
			Eigen::MatrixXcd op = j_mtxs[spin_id] * model.rho;
			double quantity = op.trace().real();
			quantities.push_back(quantity);
		}
		return quantities;
	}

	void release_observables(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);
		std::string fn;

		std::vector<double> quantities = get_quantities(model);
		for (auto spin_id = 0; spin_id < num_spins - 1; spin_id++)
		{
			fn = "j_" + std::to_string(spin_id) + model.suffix;
			save_value(quantities[spin_id], fn, save_precision);
		}
	}

	void setup_serial_data(
		Model& model,
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex) override
	{
		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);
		for (auto spin_id = 0; spin_id < num_spins - 1; spin_id++)
		{
			features_double.insert({"j_" + std::to_string(spin_id), {}});
		}
	}

	void fill_serial_data(
		Model& model,
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex) override
	{
		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);
		std::vector<double> quantities = get_quantities(model);
		for (auto spin_id = 0; spin_id < num_spins - 1; spin_id++)
		{
			features_double["j_" + std::to_string(spin_id)].push_back(quantities[spin_id]);
		}
	}
};
