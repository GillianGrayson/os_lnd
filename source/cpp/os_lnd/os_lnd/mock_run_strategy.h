#pragma once
#include "run_strategy.h"

struct MockRunStrategy : RunStrategy
{
	void run(Model& model) override
	{
		model.log_message("Mock run");
	}
};