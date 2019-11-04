#pragma once
#include "model.h"

struct ModelStrategy
{
	virtual ~ModelStrategy() = default;
	virtual void set_suffix(Model& model) = 0;
	virtual void set_sys_size(Model& model) = 0;
	virtual void set_hamiltonian(Model& model) = 0;
	virtual void set_dissipators(Model& model) = 0;
	//virtual void set_lindbladian(Model& model) = 0;
};
