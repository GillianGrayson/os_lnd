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

	void setup_aux_data(Model& model) override
	{
	}

	void setup_suffix(Model& model) override
	{
		const int name_precision = model.ini.GetInteger("global", "name_precision", 0);

		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);

		const int seed = model.ini.GetInteger("xxz", "seed", 0);

		const int diss_type = model.ini.GetInteger("xxz", "diss_type", 0);
		const auto diss_mu = model.ini.GetReal("xxz", "diss_mu", 0.0);

		const auto Delta = model.ini.GetReal("xxz", "Delta", 0.0);
		const auto h = model.ini.GetReal("xxz", "h", 0.0);

		const int quantity_index = model.ini.GetInteger("xxz", "quantity_index", 0);

		std::stringstream fns;
		fns << "_ns(" << num_spins << ")";
		fns << "_seed(" << seed << ")";
		fns << "_diss(" << diss_type << "_" << std::setprecision(name_precision) << std::fixed << diss_mu << ")";
		fns << "_prm(" << std::setprecision(name_precision) << std::fixed << Delta << "_" << h << ")";
		fns << "_q(" << quantity_index << ")";
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
		model.period = 1.0;
	}

	void setup_hamiltonian(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);

		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		
		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);
		
		const int seed = model.ini.GetInteger("xxz", "seed", 0);
		const int num_seeds = model.ini.GetInteger("xxz", "num_seeds", 0);
		
		const auto Delta = model.ini.GetReal("xxz", "Delta", 0.0);
		const auto h = model.ini.GetReal("xxz", "h", 0.0);

		model.hamiltonian = sp_mtx(model.sys_size, model.sys_size);

		energies = std::vector<double>(num_spins, 0.0);

		VSLStreamStatePtr stream;
		vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
		vslLeapfrogStream(stream, seed, num_seeds);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, num_spins, energies.data(), -h, h);

		auto fn = "energies" + model.suffix;
		save_vector(energies, fn, save_precision);

		sp_mtx sigma_0 = get_sp_eye(2);
		sp_mtx sigma_x = get_sigma_x();
		sp_mtx sigma_y = get_sigma_y();
		sp_mtx sigma_z = get_sigma_z();

		for (auto spin_id = 0; spin_id < num_spins - 1; spin_id++)
		{
			sp_mtx s_x_k0;
			sp_mtx s_x_k1 = sigma_0;
			sp_mtx s_y_k0;
			sp_mtx s_y_k1 = sigma_0;
			sp_mtx s_z_k0;
			sp_mtx s_z_k1 = sigma_0;
			if (spin_id == 0)
			{
				s_x_k0 = 0.5 * sigma_x;
				s_y_k0 = 0.5 * sigma_y;
				s_z_k0 = 0.5 * sigma_z;
			}
			else
			{
				s_x_k0 = sigma_0;
				s_y_k0 = sigma_0;
				s_z_k0 = sigma_0;
			}
			
			for (auto inner_id = 1; inner_id < num_spins; inner_id++)
			{
				if (inner_id == spin_id)
				{
					s_x_k0 = Eigen::kroneckerProduct(s_x_k0, 0.5 * sigma_x).eval();
					s_y_k0 = Eigen::kroneckerProduct(s_y_k0, 0.5 * sigma_y).eval();
					s_z_k0 = Eigen::kroneckerProduct(s_z_k0, 0.5 * sigma_z).eval();
				}
				else
				{
					s_x_k0 = Eigen::kroneckerProduct(s_x_k0, sigma_0).eval();
					s_y_k0 = Eigen::kroneckerProduct(s_y_k0, sigma_0).eval();
					s_z_k0 = Eigen::kroneckerProduct(s_z_k0, sigma_0).eval();
				}

				if (inner_id == spin_id + 1)
				{
					s_x_k1 = Eigen::kroneckerProduct(s_x_k1, 0.5 * sigma_x).eval();
					s_y_k1 = Eigen::kroneckerProduct(s_y_k1, 0.5 * sigma_y).eval();
					s_z_k1 = Eigen::kroneckerProduct(s_z_k1, 0.5 * sigma_z).eval();
				}
				else
				{
					s_x_k1 = Eigen::kroneckerProduct(s_x_k1, sigma_0).eval();
					s_y_k1 = Eigen::kroneckerProduct(s_y_k1, sigma_0).eval();
					s_z_k1 = Eigen::kroneckerProduct(s_z_k1, sigma_0).eval();
				}
			}

			model.hamiltonian += (s_x_k0 * s_x_k1 + s_y_k0 * s_y_k1 + Delta * s_z_k0 * s_z_k1 + 0.5 * energies[spin_id] * s_z_k0 + 0.5 * energies[spin_id + 1] * s_z_k1);	
		}
	}

	void setup_hamiltonian_drv(Model& model) override
	{
		model.log_message("hamiltonian_drv is absent in this model");
	}

	void setup_dissipators(Model& model) override
	{
		const int num_spins = model.ini.GetInteger("mbl", "num_spins", 0);
		const auto diss_mu = model.ini.GetReal("xxz", "diss_mu", 0.0);

		int num_diss = 4;

		sp_mtx sigma_0 = get_sp_eye(2);
		sp_mtx sigma_m = get_sigma_m();
		sp_mtx sigma_p = get_sigma_p();
		
		sp_mtx s_m_l = sigma_m;
		sp_mtx s_p_l = sigma_p;

		sp_mtx s_m_r = sigma_0;
		sp_mtx s_p_r = sigma_0;

		for (auto inner_id = 1; inner_id < num_spins; inner_id++)
		{
			s_m_l = Eigen::kroneckerProduct(s_m_l, sigma_0).eval();
			s_p_l = Eigen::kroneckerProduct(s_p_l, sigma_0).eval();

			if (inner_id == num_spins - 1)
			{
				s_m_r = Eigen::kroneckerProduct(s_m_r, sigma_m).eval();
				s_p_r = Eigen::kroneckerProduct(s_p_r, sigma_p).eval();
			}
			else
			{
				s_m_r = Eigen::kroneckerProduct(s_m_r, sigma_0).eval();
				s_p_r = Eigen::kroneckerProduct(s_p_r, sigma_0).eval();
			}
		}

		sp_mtx L1 = std::sqrt(1.0 + diss_mu) * s_p_l;
		model.dissipators.push_back(L1);
		sp_mtx L2 = std::sqrt(1.0 - diss_mu) * s_m_l;
		model.dissipators.push_back(L2);
		sp_mtx L3 = std::sqrt(1.0 - diss_mu) * s_p_r;
		model.dissipators.push_back(L3);
		sp_mtx L4 = std::sqrt(1.0 + diss_mu) * s_m_r;
		model.dissipators.push_back(L4);
	}

	void setup_lindbladian(Model& model) override
	{
		const auto diss_mu = model.ini.GetReal("xxz", "diss_mu", 0.0);

		const std::complex<double> i1(0.0, 1.0);
		const sp_mtx eye = get_sp_eye(model.sys_size);
		const sp_mtx hamiltonian_transposed(model.hamiltonian.transpose());

		model.lindbladian = -i1 * (Eigen::kroneckerProduct(eye, model.hamiltonian) - Eigen::kroneckerProduct(hamiltonian_transposed, eye));

		for (auto diss_id = 0; diss_id != model.dissipators.size(); diss_id++) 
		{
			double mult = 0.0;
			switch (diss_id)
			{
			case 0:
				mult = 1.0 + diss_mu;
				break;
			case 1:
				mult = 1.0 - diss_mu;
				break;
			case 2:
				mult = 1.0 - diss_mu;
				break;
			case 3:
				mult = 1.0 + diss_mu;
				break;
			default: 
				model.throw_error("There is only 4 dissipators in XXZ model");
			}

			sp_mtx diss = model.dissipators[diss_id];

			sp_mtx diss_tmp_1((diss.adjoint()).transpose());
			sp_mtx diss_tmp_2(diss.adjoint() * diss);
			sp_mtx diss_tmp_3(diss_tmp_2.transpose());

			model.lindbladian += 0.25 * mult * (2.0 *
				Eigen::kroneckerProduct(eye, diss) *
				Eigen::kroneckerProduct(diss_tmp_1, eye) -
				Eigen::kroneckerProduct(diss_tmp_3, eye) -
				Eigen::kroneckerProduct(eye, diss_tmp_2));
		}
	}

	void setup_lindbladian_drv(Model& model) override
	{
		model.log_message("lindbladian_drv is absent in this model");
	}

	double get_quantity(Model& model)
	{
		const int quantity_index = model.ini.GetInteger("xxz", "quantity_index", 0);
		const int num_spins = model.ini.GetInteger("xxz", "num_spins", 0);

		sp_mtx sigma_0 = get_sp_eye(2);
		sp_mtx sigma_x = get_sigma_x();
		sp_mtx sigma_y = get_sigma_y();
		sp_mtx sigma_z = get_sigma_z();

		sp_mtx s_x_k0;
		sp_mtx s_x_k1 = sigma_0;
		sp_mtx s_y_k0;
		sp_mtx s_y_k1 = sigma_0;
		sp_mtx s_z_k0;
		sp_mtx s_z_k1 = sigma_0;
		if (quantity_index == 0)
		{
			s_x_k0 = 0.5 * sigma_x;
			s_y_k0 = 0.5 * sigma_y;
			s_z_k0 = 0.5 * sigma_z;
		}
		else
		{
			s_x_k0 = sigma_0;
			s_y_k0 = sigma_0;
			s_z_k0 = sigma_0;
		}

		for (auto inner_id = 1; inner_id < num_spins; inner_id++)
		{
			if (inner_id == quantity_index)
			{
				s_x_k0 = Eigen::kroneckerProduct(s_x_k0, 0.5 * sigma_x).eval();
				s_y_k0 = Eigen::kroneckerProduct(s_y_k0, 0.5 * sigma_y).eval();
				s_z_k0 = Eigen::kroneckerProduct(s_z_k0, 0.5 * sigma_z).eval();
			}
			else
			{
				s_x_k0 = Eigen::kroneckerProduct(s_x_k0, sigma_0).eval();
				s_y_k0 = Eigen::kroneckerProduct(s_y_k0, sigma_0).eval();
				s_z_k0 = Eigen::kroneckerProduct(s_z_k0, sigma_0).eval();
			}

			if (inner_id == quantity_index + 1)
			{
				s_x_k1 = Eigen::kroneckerProduct(s_x_k1, 0.5 * sigma_x).eval();
				s_y_k1 = Eigen::kroneckerProduct(s_y_k1, 0.5 * sigma_y).eval();
				s_z_k1 = Eigen::kroneckerProduct(s_z_k1, 0.5 * sigma_z).eval();
			}
			else
			{
				s_x_k1 = Eigen::kroneckerProduct(s_x_k1, sigma_0).eval();
				s_y_k1 = Eigen::kroneckerProduct(s_y_k1, sigma_0).eval();
				s_z_k1 = Eigen::kroneckerProduct(s_z_k1, sigma_0).eval();
			}
		}

		sp_mtx j_k0 = (s_x_k0 * s_y_k1 - s_y_k0 * s_x_k1);

		Eigen::MatrixXcd op = j_k0 * model.rho;

		double quantity = model.rho.trace().real();

		return quantity;
	}

	void release_observables(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		std::string fn;

		const double quantity = get_quantity(model);
		fn = "quantity" + model.suffix;
		save_value(quantity, fn, save_precision);
	}
};
