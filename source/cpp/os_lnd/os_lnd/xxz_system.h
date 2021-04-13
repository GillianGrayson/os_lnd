#pragma once
#include "system.h"


struct XXZSystem : BaseSystem
{
	int drv_type;
	double mu;
	double ampl;
	double freq;
	double phase;
	
	XXZSystem(Model& model) : BaseSystem(model)
	{
		drv_type = model.ini.GetInteger("xxz", "drv_type", 0);
		mu = model.ini.GetReal("xxz", "mu", 0.0);
		ampl = model.ini.GetReal("xxz", "ampl", 0.0);
		freq = model.ini.GetReal("xxz", "freq", 0.0);
		phase = model.ini.GetReal("xxz", "phase", 0.0);
	}

	void operator()(const Eigen::VectorXcd& x, Eigen::VectorXcd& dxdt, const double t) override
	{
		double c0 = mu;
		double c1 = -mu;
		double c2 = -mu;
		double c3 = mu;
		
		if (drv_type == 0)
		{
			c0 = mu * (1.0 + ampl * std::sin(freq * t + phase));
			c1 = -mu * (1.0 + ampl * std::sin(freq * t + phase));
			c2 = -mu * (1.0 + ampl * std::sin(freq * t + phase));
			c3 = mu * (1.0 + ampl * std::sin(freq * t + phase));
		}
		if (drv_type == 1)
		{
			c0 = mu * (1.0 + ampl * std::sin(freq * t + phase));
			c1 = -mu * (1.0 + ampl * std::sin(freq * t + phase));
			c2 = -mu * (1.0 + ampl * std::sin(freq * t + phase));
			c3 = mu * (1.0 + ampl * std::sin(freq * t + phase));
		}

		dxdt.noalias() = (model.lindbladian + c0 * model.lindbladians_drv[0] + c1 * model.lindbladians_drv[1] + c2 * model.lindbladians_drv[2] + c3 * model.lindbladians_drv[3]) * x;
	}
};
