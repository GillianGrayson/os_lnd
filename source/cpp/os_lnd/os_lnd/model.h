#pragma once
#include <INIReader.h>
#include <vector>
#include "init.h"
#include <chrono>
#include <Eigen/Core>
#include <fmt/format.h>
#include <spdlog/spdlog.h>
#include "spdlog/sinks/stdout_color_sinks.h"
#include "save.h"

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
	std::vector<sp_mtx> dissipators;
	sp_mtx lindbladian;
	sp_mtx lindbladian_drv;
	Eigen::MatrixXcd rho;
	std::chrono::high_resolution_clock::time_point start_run_time;

	Model(INIReader& ini) : ini(ini)
	{
		start_run_time = std::chrono::high_resolution_clock::now();

		logger_type = ini.Get("global", "logger_type", "unknown");
		silent = ini.GetBoolean("global", "silent", false);
		auto logger = spdlog::stdout_color_mt(logger_type);
	}

	void throw_error(const std::string message) const
	{
		if(!silent)
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

	void log_time_duration() const
	{
		const auto run_time = std::chrono::high_resolution_clock::now();
		const auto duration = std::chrono::duration<double>(run_time - start_run_time).count();
		log_message(fmt::format("run_time = {:.16e} seconds", duration));
	}

	void log_setup_info() const
	{
		const int save_precision = ini.GetInteger("global", "save_precision", 0);
		
		log_message(fmt::format("sys_size = {}", sys_size));

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
	}
};
