#pragma once
#include "model.h"
#include "model_processor.h"

struct RunStrategy
{
	virtual ~RunStrategy() = default;
	virtual void run(
		Model& model,
		ModelProcessor& model_processor) = 0;
	virtual void run_serial(
		Model& model,
		ModelProcessor& model_processor,
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex) = 0;
};