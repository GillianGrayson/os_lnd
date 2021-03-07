#pragma once
#include "observer.h"

struct PImplObserver
{
	std::shared_ptr<BaseObserver> observer;

	PImplObserver()
	{
	}
	
	PImplObserver(std::shared_ptr<BaseObserver> obs): observer(std::move(obs))
	{
	}

	void operator()(const Eigen::VectorXcd& x, double t)
	{
		observer->operator()(x, t);
	}

	void fill_serial_features(
		std::map<std::string, std::vector<double>>& features_double,
		std::map<std::string, std::vector<std::complex<double>>>& features_complex)
	{
		observer->fill_serial_features(features_double, features_complex);
	}
};
	