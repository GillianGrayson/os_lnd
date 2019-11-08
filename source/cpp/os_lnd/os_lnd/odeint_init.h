#pragma once
#include <boost/numeric/odeint.hpp>
#include "init.h"

//typedef boost::numeric::odeint::runge_kutta_dopri5<Eigen::VectorXcd, std::complex<double>, Eigen::VectorXcd, double, boost::numeric::odeint::vector_space_algebra> error_stepper_type;
typedef boost::numeric::odeint::runge_kutta_dopri5<Eigen::VectorXcd> error_stepper_type;
typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
