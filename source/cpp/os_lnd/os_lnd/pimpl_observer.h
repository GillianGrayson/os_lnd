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
};
	