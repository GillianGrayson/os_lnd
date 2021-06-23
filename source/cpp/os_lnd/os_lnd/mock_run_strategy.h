#pragma once
#include "run_strategy.h"

struct MockRunStrategy : RunStrategy
{
	void run(Model& model, ModelProcessor& model_processor) override
	{
		model.log_message("Mock run");
	}
	
	void run_serial(
		Model& model,
		ModelProcessor& model_processor,
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex) override
	{
		model.log_message("Mock run");
	}
};