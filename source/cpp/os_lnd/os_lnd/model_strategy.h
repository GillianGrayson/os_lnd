#pragma once
#include "model.h"

struct ModelStrategy
{
	virtual ~ModelStrategy() = default;

	virtual void setup_aux_data(Model& model) = 0;
	
	virtual void setup_suffix(Model& model) = 0;
	virtual void setup_sys_size(Model& model) = 0;
	virtual void setup_period(Model& model) = 0;
	virtual void setup_hamiltonian(Model& model) = 0;
	virtual void setup_hamiltonian_drv(Model& model) = 0;
	virtual void setup_dissipators(Model& model) = 0;
	virtual void setup_lindbladian(Model& model) = 0;
	virtual void setup_lindbladians_drv(Model& model) = 0;

	virtual void release_observables(Model& model) = 0;
};
