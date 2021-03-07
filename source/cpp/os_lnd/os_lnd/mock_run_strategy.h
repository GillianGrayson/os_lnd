#pragma once
#include "run_strategy.h"

struct MockRunStrategy : RunStrategy
{
	void run(Model& model) override
	{
		model.log_message("Mock run");
	}
	
	void run_serial(
		Model& model,
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex) override
	{
		model.log_message("Mock run");
	}
};