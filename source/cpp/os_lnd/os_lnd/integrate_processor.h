#pragma once
#include "pimpl_observer.h"
#include "dimer_observer.h"
#include "mbl_observer.h"
#include "xxz_observer.h"
#include "pimpl_system.h"
#include "dimer_system.h"
#include "mbl_system.h"
#include "xxz_system.h"

struct IntegrateProcessor
{
	PImplObserver observer;
	PImplSystem system;

	Model& model;
	std::vector<double>& times;
	double& step;
	Eigen::VectorXcd& start_state;

	IntegrateProcessor(Model& model, std::vector<double>& times, double& step, Eigen::VectorXcd& start_state) : model(model), times(times), step(step), start_state(start_state)
	{
		const std::string system_type = model.ini.Get("global", "system", "unknown");
		if (system_type == "mbl")
		{
			observer = PImplObserver(std::make_shared<MBLObserver>(model, times, start_state));
			system = PImplSystem(std::make_shared<MBLSystem>(model));
		}
		else if (system_type == "dimer")
		{
			observer = PImplObserver(std::make_shared<DimerObserver>(model, times, start_state));
			system = PImplSystem(std::make_shared<DimerSystem>(model));
		}
		else if (system_type == "xxz")
		{
			observer = PImplObserver(std::make_shared<XXZObserver>(model, times, start_state));
			system = PImplSystem(std::make_shared<XXZSystem>(model));
		}
		else
		{
			model.throw_error("Unsupported system");
		}
	}

	void process() const
	{
		const runge_kutta4_stepper rk4_stepper;
		boost::numeric::odeint::integrate_times(
			rk4_stepper,
			system,
			start_state,
			times,
			step,
			observer
		);
	}

	void process_serial(
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex)
	{
		process();
		observer.fill_serial_features(features_double, features_complex);
	}
};
