#pragma once
#include "model.h"

struct ModelRunStrategy
{
	virtual ~ModelRunStrategy() = default;
	virtual void run_asymptotic_rho(Model& model) = 0;
};