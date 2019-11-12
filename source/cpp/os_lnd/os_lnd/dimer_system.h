#pragma once
#include "system.h"


struct DimerSystem : BaseSystem
{
	double E;
	int drv_type;
	double drv_ampl;
	double drv_freq;
	double drv_phase;
	
	DimerSystem(Model& model) : BaseSystem(model)
	{
		E = model.ini.GetReal("dimer", "E", 0.0);

		drv_type = model.ini.GetInteger("dimer", "drv_type", 0);
		drv_ampl = model.ini.GetReal("dimer", "drv_ampl", 0.0);
		drv_freq = model.ini.GetReal("dimer", "drv_freq", 0.0);
		drv_phase = model.ini.GetReal("dimer", "drv_phase", 0.0);
	}

	void operator()(const Eigen::VectorXcd& x, Eigen::VectorXcd& dxdt, const double t) override
	{
		double driving = 0.0;
		if (drv_type == 0)
		{
			double mod_time = std::fmod(t, model.period);
			double half_period = model.period * 0.5;
			if (mod_time < half_period)
			{
				driving = E + drv_ampl;
			}
			else
			{
				driving = E - drv_ampl;
			}
		}
		else if (drv_type == 1)
		{
			driving = E + drv_ampl * std::sin(drv_freq * t + drv_phase);
		}
		else
		{
			model.throw_error("Unsupported drv_type");
		}

		dxdt.noalias() = (model.lindbladian - driving * model.lindbladian_drv) * x;
	}
};
