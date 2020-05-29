#pragma once
#include "init.h"
#include "model.h"
#include <iomanip>
#include "save.h"
#include "model_strategy.h"
#include <unsupported/Eigen/KroneckerProduct>
#include "routines.h"

struct DimerModelStrategy : ModelStrategy
{
	void setup_aux_data(Model& model) override
	{
	}
	
	void setup_suffix(Model& model) override
	{
		const int name_precision = model.ini.GetInteger("global", "name_precision", 0);

		const int num_particles = model.ini.GetInteger("dimer", "num_particles", 0);

		const int diss_type = model.ini.GetInteger("dimer", "diss_type", 0);
		const auto diss_gamma = model.ini.GetReal("dimer", "diss_gamma", 0.0);
		
		const auto E = model.ini.GetReal("dimer", "E", 0.0);
		const auto U = model.ini.GetReal("dimer", "U", 0.0);
		const auto J = model.ini.GetReal("dimer", "J", 0.0);

		const int drv_type = model.ini.GetInteger("dimer", "drv_type", 0);
		const auto drv_ampl = model.ini.GetReal("dimer", "drv_ampl", 0.0);
		const auto drv_freq = model.ini.GetReal("dimer", "drv_freq", 0.0);
		const auto drv_phase = model.ini.GetReal("dimer", "drv_phase", 0.0);

		std::stringstream fns;
		fns << "_np(" << num_particles << ")";
		fns << "_diss(" << diss_type << "_" << std::setprecision(name_precision) << std::fixed << diss_gamma << ")";
		fns << "_prm(" << std::setprecision(name_precision) << std::fixed << E << "_" << U << "_" << J << ")";
		fns << "_drv(" << drv_type << "_" << std::setprecision(name_precision) << std::fixed << drv_ampl << "_" << drv_freq << "_" << drv_phase << ")";
		fns << ".txt";

		model.suffix = fns.str();
	}

	void setup_sys_size(Model& model) override
	{
		const int num_particles = model.ini.GetInteger("dimer", "num_particles", 0);
		
		model.sys_size = num_particles + 1;
	}

	void setup_period(Model& model) override
	{
		const auto drv_freq = model.ini.GetReal("dimer", "drv_freq", 0.0);

		const double PI = std::atan(1.0) * 4;
		
		model.period = 2.0 * PI / drv_freq;	
	}

	void setup_hamiltonian(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const int num_particles = model.ini.GetInteger("dimer", "num_particles", 0);
		auto E = model.ini.GetReal("dimer", "E", 0.0);
		auto U = model.ini.GetReal("dimer", "U", 0.0);
		U /= double(num_particles);
		const auto J = model.ini.GetReal("dimer", "J", 0.0);

		std::vector<int> rows;
		std::vector<int> cols;
		std::vector<double> vals;

		double trace = 0.0;
		for (int id = 0; id < model.sys_size; id++)
		{
			const double val_U = 2.0 * U * double(id * (id - 1) + (model.sys_size - (id + 1)) * (model.sys_size - (id + 1) - 1));
			const double val_E = E * double((model.sys_size - (id + 1)) - id);
			trace += (val_U + val_E);
			vals.push_back(val_U + val_E);
			rows.push_back(id);
			cols.push_back(id);
		}
		trace /= double(model.sys_size);

		for (int id = 0; id < model.sys_size; id++)
		{
			vals[id] -= trace;
		}

		for (int id = 0; id < (model.sys_size - 1); id++)
		{
			vals.push_back(-J * std::sqrt(double((model.sys_size - (id + 1)) * (id + 1))));
			rows.push_back(id + 1);
			cols.push_back(id);

			vals.push_back(-J * std::sqrt(double((id + 1) * (model.sys_size - (id + 1)))));
			rows.push_back(id);
			cols.push_back(id + 1);
		}
		
		std::vector<triplet> vec_triplets;
		vec_triplets.reserve(vals.size());
		for (int st_id = 0; st_id < vals.size(); st_id++)
		{
			vec_triplets.push_back(triplet(rows[st_id], cols[st_id], std::complex<double>(vals[st_id], 0.0)));
		}

		model.hamiltonian = sp_mtx(model.sys_size, model.sys_size);
		model.hamiltonian.setFromTriplets(vec_triplets.begin(), vec_triplets.end());
	}

	void setup_hamiltonian_drv(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);

		std::vector<int> rows;
		std::vector<int> cols;
		std::vector<double> vals;

		for (int id = 0; id < model.sys_size; id++)
		{
			vals.push_back(-1.0 * double((model.sys_size - (id + 1)) - id));
			rows.push_back(id);
			cols.push_back(id);
		}

		std::vector<triplet> vec_triplets;
		vec_triplets.reserve(vals.size());
		for (int st_id = 0; st_id < vals.size(); st_id++)
		{
			vec_triplets.push_back(triplet(rows[st_id], cols[st_id], std::complex<double>(vals[st_id], 0.0)));
		}

		model.hamiltonian_drv = sp_mtx(model.sys_size, model.sys_size);
		model.hamiltonian_drv.setFromTriplets(vec_triplets.begin(), vec_triplets.end());
	}

	void setup_dissipators(Model& model) override
	{
		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const int diss_type = model.ini.GetInteger("dimer", "diss_type", 0);
		
		if (diss_type == 1)
		{
			std::vector<int> rows;
			std::vector<int> cols;
			std::vector<double> vals;

			for (int id = 0; id < model.sys_size; id++)
			{
				vals.push_back(double(id - (model.sys_size - (id + 1))));
				rows.push_back(id);
				cols.push_back(id);
			}

			for (int id = 0; id < (model.sys_size - 1); id++)
			{
				vals.push_back(-std::sqrt(double((model.sys_size - (id + 1)) * (id + 1))));
				rows.push_back(id + 1);
				cols.push_back(id);

				vals.push_back(+std::sqrt(double((id + 1) * (model.sys_size - (id + 1)))));
				rows.push_back(id);
				cols.push_back(id + 1);
			}

			std::vector<triplet> vec_triplets;
			vec_triplets.reserve(vals.size());
			for (int id = 0; id < vals.size(); id++)
			{
				vec_triplets.push_back(triplet(rows[id], cols[id], std::complex<double>(vals[id], 0.0)));
			}

			sp_mtx mtx(model.sys_size, model.sys_size);
			mtx.setFromTriplets(vec_triplets.begin(), vec_triplets.end());

			model.dissipators.push_back(mtx);
		}
		else
		{
			model.throw_error("Unsupported dissipator type");
		}
	}

	void setup_lindbladian(Model& model) override
	{
		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const int num_particles = model.ini.GetInteger("dimer", "num_particles", 0);
		double diss_gamma = model.ini.GetReal("dimer", "diss_gamma", 0.0);
		diss_gamma /= double(num_particles);

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

	void setup_lindbladian_drv(Model& model) override
	{
		const auto debug_dump = model.ini.GetBoolean("global", "debug_dump", false);
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);

		const std::complex<double> i1(0.0, 1.0);
		const sp_mtx eye = get_sp_eye(model.sys_size);
		const sp_mtx hamiltonian_transposed(model.hamiltonian_drv.transpose());

		model.lindbladian_drv = -i1 * (Eigen::kroneckerProduct(eye, model.hamiltonian_drv) - Eigen::kroneckerProduct(hamiltonian_transposed, eye));
	}

	void release_observables(Model& model) override
	{
	}
};
