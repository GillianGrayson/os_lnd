#pragma once
#include "system.h"


struct XXZSystem : BaseSystem
{
	int drv_type;
	double mu;
	double T1;
	double T2;
	
	XXZSystem(Model& model) : BaseSystem(model)
	{
		drv_type = model.ini.GetInteger("xxz", "drv_type", 0);
		mu = model.ini.GetReal("xxz", "mu", 0.0);
		T1 = model.ini.GetReal("xxz", "T1", 0.0);
		T2 = model.ini.GetReal("xxz", "T2", 0.0);
	}

	void operator()(const Eigen::VectorXcd& x, Eigen::VectorXcd& dxdt, const double t) override
	{
		double c0 = mu;
		double c1 = -mu;
		double c2 = -mu;
		double c3 = mu;
		
		if (drv_type == 1)
		{
			c0 = mu * std::cos(std::fmod(t, T1) / T1);
			c1 = -mu * std::cos(std::fmod(t, T1) / T1);
			c2 = -mu * std::cos(std::fmod(t, T2) / T2);
			c3 = mu * std::cos(std::fmod(t, T2) / T2);
		}

		dxdt.noalias() = (model.lindbladian + c0 * model.lindbladians_drv[0] + c1 * model.lindbladians_drv[1] + c2 * model.lindbladians_drv[2] + c3 * model.lindbladians_drv[3]) * x;
	}
};
