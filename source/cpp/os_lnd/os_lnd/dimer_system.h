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
			double sinus = std::sin(drv_freq * t + drv_phase);
			if (sinus > 0)
			{
				driving = drv_ampl;
			}
			else
			{
				driving = -drv_ampl;
			}
		}
		else if (drv_type == 1)
		{
			driving = drv_ampl * std::sin(drv_freq * t + drv_phase);
		}
		else
		{
			model.throw_error("Unsupported drv_type");
		}

		dxdt.noalias() = (model.lindbladian - driving * model.lindbladians_drv[0]) * x;
	}
};
