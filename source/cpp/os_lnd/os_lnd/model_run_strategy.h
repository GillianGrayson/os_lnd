#pragma once
#include "model.h"

struct ModelRunStrategy
{
	virtual ~ModelRunStrategy() = default;
	virtual void run(Model& model) = 0;
};