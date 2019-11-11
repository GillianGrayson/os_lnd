#pragma once
#include <INIReader.h>
#include <vector>
#include "init.h"
#include <chrono>
#include <Eigen/Core>
#include <fmt/format.h>
#include <spdlog/spdlog.h>
#include "spdlog/sinks/stdout_color_sinks.h"

struct Model
{
	INIReader ini;
	std::string logger_type;
	bool silent;

	std::string suffix;
	int sys_size;
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

	void throw_error(const std::string message)
	{
		if(!silent)
		{
			spdlog::get(logger_type)->error(message);
		}
		throw std::runtime_error(message);
	}

	void log_message(const std::string message)
	{
		if (!silent)
		{
			spdlog::get(logger_type)->info(message);
		}
	}

	void log_time_duration()
	{
		const auto curr_time = std::chrono::high_resolution_clock::now();
		const auto duration = std::chrono::duration<double>(curr_time - start_run_time).count();
		log_message(fmt::format("run_time = {:.16e} seconds", duration));
	}
};
