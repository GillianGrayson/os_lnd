#pragma once
#include <boost/numeric/odeint.hpp>
#include "init.h"

using namespace boost::numeric::odeint;

typedef runge_kutta_dopri5<Eigen::VectorXcd> error_stepper_type;
typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;