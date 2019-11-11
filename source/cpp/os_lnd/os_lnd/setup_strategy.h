#pragma once
#include "model.h"

struct SetupStrategy
{
	virtual ~SetupStrategy() = default;
	virtual void setup_suffix(Model& model) = 0;
	virtual void setup_sys_size(Model& model) = 0;
	virtual void setup_hamiltonian(Model& model) = 0;
	virtual void setup_hamiltonian_drv(Model& model) = 0;
	virtual void setup_dissipators(Model& model) = 0;
	virtual void setup_lindbladian(Model& model) = 0;
	virtual void setup_lindbladian_drv(Model& model) = 0;
};
