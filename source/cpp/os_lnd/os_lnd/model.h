#pragma once
#include <INIReader.h>
#include <vector>
#include "init.h"
#include <chrono>
#include <Eigen/Core>
#include <fmt/format.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include "save.h"
#include "memory_usage.h"
#include "routines.h"

struct Model
{
	INIReader ini;
	std::string logger_type;
	bool silent;

	std::string suffix;
	int sys_size;
	double period;
	sp_mtx hamiltonian;
	sp_mtx hamiltonian_drv;
	ds_mtx hamiltonian_dense;
	ds_mtx hamiltonian_drv_dense;
	std::vector<sp_mtx> dissipators;
	std::vector<ds_mtx> dissipators_dense;
	sp_mtx lindbladian;
	sp_mtx lindbladian_drv;
	ds_mtx lindbladian_dense;
	ds_mtx lindbladian_drv_dense;
	std::vector<sp_mtx> f_basis;
	Eigen::MatrixXcd rho;
	std::chrono::high_resolution_clock::time_point start_run_time;
	std::vector<double> run_times;

	Model(INIReader& ini) : ini(ini)
	{
		start_run_time = std::chrono::high_resolution_clock::now();

		logger_type = ini.Get("global", "logger_type", "unknown");
		silent = ini.GetBoolean("global", "silent", false);
		auto logger = spdlog::stdout_color_mt(logger_type);
	}

	void throw_error(const std::string message) const
	{
		if (!silent)
		{
			spdlog::get(logger_type)->error(message);
		}
		throw std::runtime_error(message);
	}

	void log_message(const std::string message) const
	{
		if (!silent)
		{
			spdlog::get(logger_type)->info(message);
		}
	}

	void log_time_duration()
	{
		const int save_precision = ini.GetInteger("global", "save_precision", 0);
		
		const auto run_time = std::chrono::high_resolution_clock::now();
		double duration = std::chrono::duration<double>(run_time - start_run_time).count();
		log_message(fmt::format("run_time = {:.16e} seconds", duration));

		run_times.push_back(duration);

		auto fn = "run_times" + suffix;
		save_vector(run_times, fn, save_precision);
	}

	void log_setup_info()
	{
		const int save_precision = ini.GetInteger("global", "save_precision", 0);
		
		log_message(fmt::format("sys_size = {}\n", sys_size));

		const auto lindbladian_non_zeros = lindbladian.nonZeros();
		const auto lindbladian_non_zeros_part = double(lindbladian_non_zeros) / (std::pow(double(sys_size), 4.0));
		if (lindbladian_non_zeros > 0)
		{
			log_message(fmt::format("Number of non-zero elements in lindbladian = {}", lindbladian_non_zeros));
			log_message(fmt::format("Part of non-zero elements in lindbladian = {:.16e}\n", lindbladian_non_zeros_part));
		}

		const auto lindbladian_drv_non_zeros = lindbladian_drv.nonZeros();
		const auto lindbladian_drv_non_zeros_part = double(lindbladian_drv_non_zeros) / (std::pow(double(sys_size), 4.0));
		if (lindbladian_drv_non_zeros > 0)
		{
			log_message(fmt::format("Number of non-zero elements in lindbladian_drv = {}", lindbladian_drv_non_zeros));
			log_message(fmt::format("Part of non-zero elements in lindbladian_drv = {:.16e}\n", lindbladian_drv_non_zeros_part));
		}

		std::vector<double> non_zeros_parts;
		non_zeros_parts.push_back(lindbladian_non_zeros_part);
		non_zeros_parts.push_back(lindbladian_drv_non_zeros_part);

		auto fn = "non_zeros_parts" + suffix;
		save_vector(non_zeros_parts, fn, save_precision);

		log_time_duration();
		log_memory_usage();
	}

	void log_memory_usage() const
	{
		const int save_precision = ini.GetInteger("global", "save_precision", 0);

		double currentSize = double(getCurrentRSS()) / std::pow(1024.0, 2);
		double peakSize = double(getPeakRSS()) / std::pow(1024.0, 2);

		log_message(fmt::format("Current RSS (physical memory use) = {:.2f} Mb", currentSize));
		log_message(fmt::format("Peak RSS (physical memory use) = {:.2f} Mb\n", peakSize));

		std::vector<double> mem_info;
		mem_info.push_back(currentSize);
		mem_info.push_back(peakSize);

		auto fn = "mem_info" + suffix;
		save_vector(mem_info, fn, save_precision);
	}

	void save_data() const
	{
		const int save_precision = ini.GetInteger("global", "save_precision", 0);
		const auto debug_dump = ini.GetBoolean("global", "debug_dump", false);
		const auto save_hamiltonians = ini.GetBoolean("global", "save_hamiltonians", false);
		const auto save_dissipators = ini.GetBoolean("global", "save_dissipators", false);
		const auto save_lindbladians = ini.GetBoolean("global", "save_lindbladians", false);
		const auto save_f_basis = ini.GetBoolean("global", "save_f_basis", false);

		if (debug_dump || save_hamiltonians)
		{
			auto is_sp_empty = (hamiltonian.outerSize() > 0) ? false : true;
			auto is_ds_empty = (hamiltonian_dense.outerSize() > 0) ? false : true;
			auto fn = "hamiltonian_mtx" + suffix;
			if (!is_sp_empty)
			{
				save_sp_mtx(hamiltonian, fn, save_precision);
			}
			else if (is_sp_empty and !is_ds_empty)
			{
				save_dense_mtx(hamiltonian_dense, fn, save_precision);
			}
			else
			{
				log_message("Can't save: Empty hamiltonian");
			}

			is_sp_empty = (hamiltonian_drv.outerSize() > 0) ? false : true;
			is_ds_empty = (hamiltonian_drv_dense.outerSize() > 0) ? false : true;
			fn = "hamiltonian_drv_mtx" + suffix;
			if (!is_sp_empty)
			{
				save_sp_mtx(hamiltonian_drv, fn, save_precision);
			}
			else if (is_sp_empty and !is_ds_empty)
			{
				save_dense_mtx(hamiltonian_drv_dense, fn, save_precision);
			}
			else
			{
				log_message("Can't save: Empty hamiltonian_drv");
			}
		}

		if (debug_dump || save_dissipators)
		{
			for (auto diss_id = 0; diss_id < dissipators.size(); diss_id++)
			{
				auto is_sp_empty = (dissipators[diss_id].outerSize() > 0) ? false : true;
				auto is_ds_empty = (dissipators_dense[diss_id].outerSize() > 0) ? false : true;
				auto fn = fmt::format("diss_{:d}_mtx", diss_id) + suffix;
				if (!is_sp_empty)
				{
					save_sp_mtx(dissipators[diss_id], fn, save_precision);
				}
				else if (is_sp_empty and !is_ds_empty)
				{
					save_dense_mtx(dissipators_dense[diss_id], fn, save_precision);
				}
				else
				{
					log_message("Can't save: Empty dissipators");
				}
			}
		}

		if (debug_dump || save_lindbladians)
		{
			auto is_sp_empty = (lindbladian.outerSize() > 0) ? false : true;
			auto is_ds_empty = (lindbladian_dense.outerSize() > 0) ? false : true;
			auto fn = "lindbladian_mtx" + suffix;

			if (!is_sp_empty)
			{
				save_sp_mtx(lindbladian, fn, save_precision);
			}
			else if (is_sp_empty and !is_ds_empty)
			{
				save_dense_mtx(lindbladian_dense, fn, save_precision);
			}
			else
			{
				log_message("Can't save: Empty lindbladian");
			}
			
			is_sp_empty = (lindbladian_drv.outerSize() > 0) ? false : true;
			is_ds_empty = (lindbladian_drv_dense.outerSize() > 0) ? false : true;
			fn = "lindbladian_drv_mtx" + suffix;

			if (!is_sp_empty)
			{
				save_sp_mtx(lindbladian_drv, fn, save_precision);
			}
			else if (is_sp_empty and !is_ds_empty)
			{
				save_dense_mtx(lindbladian_drv_dense, fn, save_precision);
			}
			else
			{
				log_message("Can't save: Empty lindbladian_drv");
			}
		}

		if (debug_dump || save_f_basis)
		{
			for (auto fb_id = 0; fb_id < f_basis.size(); fb_id++)
			{
				auto fn = fmt::format("f_basis_{:d}_mtx", fb_id) + suffix;
				save_sp_mtx(f_basis[fb_id], fn, save_precision);
			}
		}
	}

	void init_f_basis()
	{
		const int n = sys_size;

		sp_mtx tmp = get_sp_eye(n) / std::sqrt(double(n));
		f_basis.push_back(tmp);

		log_message("F-basis init...");
		log_time_duration();

		double sqrt_2 = std::sqrt(2.0);

		for (auto i = 0; i < n; i++)
		{
			for (auto j = i + 1; j < n; j++)
			{
				std::vector<triplet> vec_triplets(2);
				vec_triplets[0] = triplet(i, j, std::complex<double>(1.0 / sqrt_2, 0.0));
				vec_triplets[1] = triplet(j, i, std::complex<double>(1.0 / sqrt_2, 0.0));
				tmp = sp_mtx(n, n);
				tmp.setFromTriplets(vec_triplets.begin(), vec_triplets.end());
				f_basis.push_back(tmp);

				vec_triplets[0] = triplet(i, j, std::complex<double>(0.0, -1.0 / sqrt_2));
				vec_triplets[1] = triplet(j, i, std::complex<double>(0.0, 1.0 / sqrt_2));
				tmp = sp_mtx(n, n);
				tmp.setFromTriplets(vec_triplets.begin(), vec_triplets.end());
				f_basis.push_back(tmp);
			}
		}

		for (auto i = 0; i < n - 1; i++)
		{
			std::vector<triplet> vec_triplets(i + 2);
			for (auto j = 0; j < i + 1; j++)
			{
				vec_triplets[j] = triplet(j, j, std::complex<double>(1.0 / std::sqrt(double((i + 1) * (i + 2))), 0.0));
			}
			vec_triplets[i + 1] = triplet(i + 1, i + 1, std::complex<double>(-double(i + 1) / std::sqrt(double((i + 1) * (i + 2))), 0.0));
			tmp = sp_mtx(n, n);
			tmp.setFromTriplets(vec_triplets.begin(), vec_triplets.end());
			f_basis.push_back(tmp);
		}
		
		log_message("F-basis init complete");
		log_time_duration();
	}
};
