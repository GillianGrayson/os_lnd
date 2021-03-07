#pragma once
#include "model.h"

struct RunStrategy
{
	virtual ~RunStrategy() = default;
	virtual void run(Model& model) = 0;
	virtual void run_serial(
		Model& model,
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex) = 0;
};