#pragma once
#include "model.h"

struct RunStrategy
{
	virtual ~RunStrategy() = default;
	virtual void run(Model& model) = 0;
};