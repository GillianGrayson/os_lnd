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


inline int bit_count(int value)
{
	auto count = 0;

	while (value > 0)
	{
		if ((value & 1) == 1)
		{
			count++;
		}
		value >>= 1;
	}

	return count;
}

inline int bit_at(const int value, const int position)
{
	if ((value >> position) & 1)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

inline std::vector<int> convert_int_to_vector_of_bits(int x, const int size)
{
	std::vector<int> res;

	auto id = 0;

	while (id < size)
	{
		if (x & 1)
		{
			res.push_back(1);
		}
		else
		{
			res.push_back(0);
		}

		x >>= 1;
		id++;
	}

	reverse(res.begin(), res.end());

	return res;
}

struct MBLModelStrategy : ModelStrategy
{
	std::vector<int> adj;
	std::vector<int> x_to_id;
	std::vector<int> id_to_x;
	std::vector<double> energies;

	void setup_aux_data(Model& model) override
	{
		const int num_spins = model.ini.GetInteger("mbl", "num_spins", 0);
		const int num_global_states = int(std::pow(2, num_spins));
		adj = std::vector<int>(num_global_states, 0);
		x_to_id = std::vector<int>(num_global_states, 0);
		id_to_x = std::vector<int>(model.sys_size, 0);

		int state_id = 0;
		for (int g_state_id = 0; g_state_id < num_global_states; g_state_id++)
		{
			if ((bit_count(g_state_id) == 2) && (bit_count(g_state_id & (g_state_id << 1)) == 1))
			{
				adj[g_state_id] = 1;
			}
			else
			{
				adj[g_state_id] = 0;
			}

			if (bit_count(g_state_id) == num_spins / 2)
			{
				x_to_id[g_state_id] = state_id + 1;
				id_to_x[state_id] = g_state_id;
				state_id++;
			}
			else
			{
				x_to_id[g_state_id] = 0;
			}
		}

		if (state_id != model.sys_size)
		{
			model.throw_error("Something wrong with MBL generation");
		}
	}
	
	void setup_suffix(Model& model) override
	{
		const int name_precision = model.ini.GetInteger("global", "name_precision", 0);

		const int num_spins = model.ini.GetInteger("mbl", "num_spins", 0);

		const int seed = model.ini.GetInteger("mbl", "seed", 0);

		const int diss_type = model.ini.GetInteger("mbl", "diss_type", 0);
		const auto diss_phase = model.ini.GetReal("mbl", "diss_phase", 0.0);
		const auto diss_gamma = model.ini.GetReal("mbl", "diss_gamma", 0.0);

		const auto W = model.ini.GetReal("mbl", "W", 0.0);
		const auto U = model.ini.GetReal("mbl", "U", 0.0);
		const auto J = model.ini.GetReal("mbl", "J", 0.0);

		std::stringstream fns;
		fns << "_ns(" << num_spins << ")";
		fns << "_seed(" << seed << ")";
		fns << "_diss(" << diss_type << "_" << std::setprecision(name_precision) << std::fixed << diss_phase << "_" << diss_gamma << ")";
		fns << "_prm(" << std::setprecision(name_precision) << std::fixed << W << "_" << U << "_" << J << ")";
		fns << ".txt";

		model.suffix = fns.str();
	}

	void setup_sys_size(Model& model) override
	{
		const int num_spins = model.ini.GetInteger("mbl", "num_spins", 0);
		if (std::div(num_spins, 2).rem != 0)
		{
			model.throw_error("num_spins must be divided by 2 without remainder");
		}
		
		model.sys_size = gcem::binomial_coef(num_spins, num_spins / 2);
	}

	void setup_period(Model& model) override
	{
		model.period = 1.0;
	}

	void setup_hamiltonian(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);

		model.hamiltonian = get_disorder_mtx(model) + get_interaction_mtx(model) + get_hopping_mtx(model);
	}

	void setup_hamiltonian_drv(Model& model) override
	{
		model.log_message("hamiltonian_drv is absent in this model");
	}

	void setup_dissipators(Model& model) override
	{
		const int num_spins = model.ini.GetInteger("mbl", "num_spins", 0);
		const int diss_type = model.ini.GetInteger("mbl", "diss_type", 0);
		int num_diss;
		if (diss_type == 0)
		{
			num_diss = num_spins;
		}
		else if (diss_type == 1)
		{
			num_diss = num_spins - 1;
		}
		else
		{
			model.throw_error("Unsupported dissipator type");
		}

		for (int diss_id = 0; diss_id < num_diss; diss_id++)
		{
			sp_mtx diss;
			if (diss_type == 0)
			{
				diss = get_diss_mtx_type_0(model, diss_id);
			}
			else if (diss_type == 1)
			{
				diss = get_diss_mtx_type_1(model, diss_id);
			}
			else
			{
				model.throw_error("Unsupported dissipator type");
			}

			model.dissipators.push_back(diss);	
		}
	}

	void setup_lindbladian(Model& model) override
	{
		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const auto diss_gamma = model.ini.GetReal("mbl", "diss_gamma", 0.0);

		const std::complex<double> i1(0.0, 1.0);
		const sp_mtx eye = get_sp_eye(model.sys_size);
		const sp_mtx hamiltonian_transposed(model.hamiltonian.transpose());
		
		model.lindbladian = -i1 * (Eigen::kroneckerProduct(eye, model.hamiltonian) - Eigen::kroneckerProduct(hamiltonian_transposed, eye));

		for (const auto& diss : model.dissipators)
		{
			sp_mtx diss_tmp_1((diss.adjoint()).transpose());
			sp_mtx diss_tmp_2(diss.adjoint() * diss);
			sp_mtx diss_tmp_3(diss_tmp_2.transpose());

			model.lindbladian += 0.5 * diss_gamma * (2.0 * 
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

	sp_mtx get_disorder_mtx(Model& model)
	{
		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const int num_spins = model.ini.GetInteger("mbl", "num_spins", 0);
		const int seed = model.ini.GetInteger("mbl", "seed", 0);
		const int num_seeds = model.ini.GetInteger("mbl", "num_seeds", 0);
		const auto W = model.ini.GetReal("mbl", "W", 0.0);

		energies = std::vector<double>(num_spins, 0.0);

		VSLStreamStatePtr stream;
		vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
		vslLeapfrogStream(stream, seed, num_seeds);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, num_spins, energies.data(), -1.0, 1.0);

		auto fn = "energies" + model.suffix;
		save_vector(energies, fn, save_precision);

		std::vector<triplet> vec_triplets;
		vec_triplets.reserve(model.sys_size);
		for (int st_id = 0; st_id < model.sys_size; st_id++)
		{
			std::vector<int> vb = convert_int_to_vector_of_bits(id_to_x[st_id], num_spins);
			double sum = 0.0;
			for (int cell_id = 0; cell_id < num_spins; cell_id++)
			{
				sum += double(vb[cell_id]) * energies[cell_id];
			}
			sum *= 1.0 * W;

			vec_triplets.push_back(triplet(st_id, st_id, std::complex<double>(sum, 0.0)));
		}

		sp_mtx mtx(model.sys_size, model.sys_size);
		mtx.setFromTriplets(vec_triplets.begin(), vec_triplets.end());

		if (debug_dump)
		{
			fn = "disorder_mtx" + model.suffix;
			save_sp_mtx(mtx, fn, save_precision);
		}
		
		return mtx;
	}

	sp_mtx get_interaction_mtx(Model& model)
	{
		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const auto U = model.ini.GetReal("mbl", "U", 0.0);

		std::vector<triplet> vec_triplets;
		vec_triplets.reserve(model.sys_size);
		for (int st_id = 0; st_id < model.sys_size; st_id++)
		{
			double val = U * bit_count(id_to_x[st_id] & (id_to_x[st_id] << 1));

			vec_triplets.push_back(triplet(st_id, st_id, std::complex<double>(val, 0.0)));
		}

		sp_mtx mtx(model.sys_size, model.sys_size);
		mtx.setFromTriplets(vec_triplets.begin(), vec_triplets.end());

		if (debug_dump)
		{
			auto fn = "interaction_mtx" + model.suffix;
			save_sp_mtx(mtx, fn, save_precision);
		}

		return mtx;
	}

	sp_mtx get_hopping_mtx(Model& model)
	{
		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const auto J = model.ini.GetReal("mbl", "J", 0.0);

		std::vector<int> rows;
		std::vector<int> cols;
		std::vector<double> vals;

		for (int id_1 = 0; id_1 < model.sys_size; id_1++)
		{
			for (int id_2 = 0; id_2 < model.sys_size; id_2++)
			{
				if (adj[id_to_x[id_1] ^ id_to_x[id_2]] > 0)
				{
					vals.push_back(-J);
					rows.push_back(id_1);
					cols.push_back(id_2);
				}
			}
		}

		std::vector<triplet> vec_triplets;
		vec_triplets.reserve(vals.size());
		for (int st_id = 0; st_id < vals.size(); st_id++)
		{
			vec_triplets.push_back(triplet(rows[st_id], cols[st_id], std::complex<double>(vals[st_id], 0.0)));
		}

		sp_mtx mtx(model.sys_size, model.sys_size);
		mtx.setFromTriplets(vec_triplets.begin(), vec_triplets.end());

		if (debug_dump)
		{
			auto fn = "hopping_mtx" + model.suffix;
			save_sp_mtx(mtx, fn, save_precision);
		}

		return mtx;
	}

	sp_mtx get_diss_mtx_type_0(Model& model, const int diss_id)
	{
		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);

		std::vector<triplet> vec_triplets;
		vec_triplets.reserve(model.sys_size);
		for (int st_id = 0; st_id < model.sys_size; st_id++)
		{
			double val = double(bit_at(id_to_x[st_id], diss_id));

			vec_triplets.push_back(triplet(st_id, st_id, std::complex<double>(val, 0.0)));
		}

		sp_mtx mtx(model.sys_size, model.sys_size);
		mtx.setFromTriplets(vec_triplets.begin(), vec_triplets.end());

		return mtx;
	}

	sp_mtx get_diss_mtx_type_1(Model& model, const int diss_id)
	{
		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const int num_spins = model.ini.GetInteger("mbl", "num_spins", 0);
		const auto diss_phase = model.ini.GetReal("mbl", "diss_phase", 0.0);

		std::vector<int> rows;
		std::vector<int> cols;
		std::vector<std::complex<double>> vals;

		for (int st_id_1 = 0; st_id_1 < model.sys_size; st_id_1++)
		{
			std::complex<double> val(double(bit_at(id_to_x[st_id_1], diss_id)) - double(bit_at(id_to_x[st_id_1], diss_id + 1)), 0.0);
			vals.push_back(val);
			rows.push_back(st_id_1);
			cols.push_back(st_id_1);
		}

		int tmp;
		for (int state_id_1 = 0; state_id_1 < model.sys_size; state_id_1++)
		{
			int row_id = state_id_1;

			tmp = bit_at(id_to_x[state_id_1], diss_id) - bit_at(id_to_x[state_id_1], diss_id + 1);

			int col_id = 0;

			if (tmp == 0)
			{
				col_id = state_id_1;
			}
			else
			{
				for (int state_id_2 = 0; state_id_2 < model.sys_size; state_id_2++)
				{
					if (adj[id_to_x[state_id_1] ^ id_to_x[state_id_2]])
					{
						std::vector<int> adjacency_bits = convert_int_to_vector_of_bits(id_to_x[state_id_1] ^ id_to_x[state_id_2], num_spins);
						std::vector<int> hop;
						for (int cell_id = 0; cell_id < num_spins; cell_id++)
						{
							if (adjacency_bits[cell_id])
							{
								hop.push_back(cell_id);
							}
						}

						for (int ad_cell_id = 0; ad_cell_id < hop.size(); ad_cell_id++)
						{
							hop[ad_cell_id] = (num_spins - 1) - hop[ad_cell_id];
						}

						if (hop[1] == diss_id)
						{
							if (bit_at(id_to_x[state_id_1], diss_id))
							{
								col_id = state_id_2;
								std::complex<double> val(-cos(-diss_phase), -sin(-diss_phase));

								vals.push_back(val);
								rows.push_back(row_id);
								cols.push_back(col_id);
							}
							else
							{
								col_id = state_id_2;
								std::complex<double> val(cos(diss_phase), sin(diss_phase));
								vals.push_back(val);
								rows.push_back(row_id);
								cols.push_back(col_id);
							}
						}

						adjacency_bits.clear();
						hop.clear();
					}
				}
			}
		}

		std::vector<triplet> vec_triplets;
		vec_triplets.reserve(vals.size());
		for (int st_id = 0; st_id < vals.size(); st_id++)
		{
			vec_triplets.push_back(triplet(rows[st_id], cols[st_id], vals[st_id]));
		}

		sp_mtx mtx(model.sys_size, model.sys_size);
		mtx.setFromTriplets(vec_triplets.begin(), vec_triplets.end());

		return mtx;
	}

	double get_ratio(Model& model)
	{
		Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
		ces.compute(model.rho);
		auto evals_raw = ces.eigenvalues();
		std::vector<double> evals(model.sys_size);
		for (int st_id = 0; st_id < model.sys_size; st_id ++)
		{
			evals[st_id] = evals_raw[st_id].real();
		}
		std::sort(evals.begin(), evals.end());

		int skip = model.sys_size / 3;

		std::vector<double> cutted_evals = std::vector<double>(evals.begin() + skip, evals.end() - skip);
		int size_cutted = int(cutted_evals.size());
		
		std::vector<double> evals_diff(size_cutted - 1);
		for (int st_id = 0; st_id < size_cutted - 1; st_id++)
		{
			evals_diff[st_id] = cutted_evals[st_id + 1] - cutted_evals[st_id];
			if (evals_diff[st_id] < std::numeric_limits<double>::epsilon())
			{
				model.log_message("Can't calculate ratio\n");
				return 0.0;
			}
		}
		
		std::vector<double> evals_r(size_cutted - 2);
		for (int st_id = 0; st_id < size_cutted - 2; st_id++)
		{
			evals_r[st_id] = evals_diff[st_id + 1] / evals_diff[st_id];
			double straight = evals_r[st_id];
			double inverse = 1.0 / evals_r[st_id];
			evals_r[st_id] = std::min(straight, inverse);	
		}

		double sum = std::accumulate(evals_r.begin(), evals_r.end(), 0.0);
		double ratio = sum / evals_r.size();

		return ratio;
	}

	double get_entanglement_entropy(Model& model)
	{
		const int num_spins = model.ini.GetInteger("mbl", "num_spins", 0);
		int xtd_size = int(std::pow(2, num_spins));
		
		Eigen::MatrixXcd rho_xtd = Eigen::MatrixXcd::Zero(xtd_size, xtd_size);
		for (int st_id_1 = 0; st_id_1 < model.sys_size; st_id_1++)
		{
			for (int st_id_2 = 0; st_id_2 < model.sys_size; st_id_2++)
			{
				int xtd_st_id_1 = id_to_x[st_id_1];
				int xtd_st_id_2 = id_to_x[st_id_2];

				std::complex<double> val(model.rho(st_id_1, st_id_2).real(), model.rho(st_id_1, st_id_2).imag());
				rho_xtd(xtd_st_id_1, xtd_st_id_2) = val;
			}
		}

		int red_size = int(std::pow(2, num_spins / 2));
		double eps_eval = 1.0e-14;

		Eigen::MatrixXcd R = Eigen::MatrixXcd::Zero(xtd_size, xtd_size);

		for (int k1 = 0; k1 < red_size; k1++)
		{
			for (int k2 = 0; k2 < red_size; k2++)
			{
				for (int s1 = 0; s1 < red_size; s1++)
				{
					for (int s2 = 0; s2 < red_size; s2++)
					{
						int R_id_1 = k1 * red_size + k2;
						int R_id_2 = s1 * red_size + s2;

						int rho_id_1 = k1 * red_size + s1;
						int rho_id_2 = k2 * red_size + s2;

						std::complex<double> val(rho_xtd(rho_id_1, rho_id_2).real(), rho_xtd(rho_id_1, rho_id_2).imag());
						R(R_id_1, R_id_2) = val;
					}
				}
			}
		}

		Eigen::JacobiSVD<Eigen::MatrixXcd> svd(R);
		auto svals = svd.singularValues();

		double sum = 0.0;
		for (int xtd_st_id = 0; xtd_st_id < xtd_size; xtd_st_id++)
		{
			sum += svals[xtd_st_id] * svals[xtd_st_id];
		}

		double ee = 0.0; // entanglement entropy
		for (int xtd_st_id = 0; xtd_st_id < xtd_size; xtd_st_id++)
		{
			const auto curr_mu = svals[xtd_st_id] * svals[xtd_st_id] / sum;
			if (fabs(curr_mu) > eps_eval)
			{
				ee -= curr_mu * log2(curr_mu);
			}
		}
		
		return ee;
	}

	double get_imbalance(Model& model)
	{
		const int num_spins = model.ini.GetInteger("mbl", "num_spins", 0);
		
		std::vector<double> n_part(num_spins);

		for (int st_id = 0; st_id < model.sys_size; st_id++)
		{
			std::vector<int> vb = convert_int_to_vector_of_bits(id_to_x[st_id], num_spins);

			for (int cell_id = 0; cell_id < num_spins; cell_id++)
			{
				n_part[cell_id] += model.rho(st_id, st_id).real() * double(vb[cell_id]);
			}
		}

		double sum_odd = 0.0;
		double sum_even = 0.0;
		double sum_all = 0.0;

		for (int cell_id = 0; cell_id < num_spins; cell_id++)
		{
			if (cell_id % 2 == 0)
			{
				sum_odd += n_part[cell_id];
			}
			else
			{
				sum_even += n_part[cell_id];
			}
		}
		sum_all = sum_even + sum_odd;

		double imbalance = (sum_odd - sum_even) / sum_all;

		return imbalance;
	}

	void release_observables(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		std::string fn;

		const double ratio = get_ratio(model);
		fn = "ratio" + model.suffix;
		save_value(ratio, fn, save_precision);

		const double ee = get_entanglement_entropy(model);
		fn = "ee" + model.suffix;
		save_value(ee, fn, save_precision);

		const double imbalance = get_imbalance(model);
		fn = "imbalance" + model.suffix;
		save_value(imbalance, fn, save_precision);
	}

	void setup_serial_data(
		Model& model,
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex) override
	{
	}

	void fill_serial_data(
		Model& model,
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex) override
	{
	}
};
