//#include "model_setup_processor.h"
//#include "model_run_processor.h"
//
//
//int main(int argc, char* argv[])
//{
//	INIReader ini("config.ini");
//	if (ini.ParseError() < 0)
//	{
//		throw std::runtime_error("Can't load 'test.ini'\n");
//	}
//
//	Model model(ini);
//
//	ModelSetupProcessor setup_processor;
//	setup_processor.set_model_strategy(model);
//	setup_processor.process_model(model);
//
//	ModelRunProcessor run_processor;
//	run_processor.set_model_strategy(model);
//	run_processor.process_model(model);
//}

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include <Eigen/Core> 
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

// define vector_space_algebra for Eigen::Matrix
namespace boost::numeric::odeint {
	template<typename B, int S1, int S2, int O, int M1, int M2>
	struct algebra_dispatcher< Eigen::Matrix<B, S1, S2, O, M1, M2> > {
		typedef vector_space_algebra algebra_type;
	};
}

// define abs() for Eigen::Matrix
namespace Eigen {
	template<typename D, int Rows, int Cols>
	Matrix<D, Rows, Cols> abs(Matrix<D, Rows, Cols> const& m) {
		return m.cwiseAbs();
	}
}

typedef Eigen::Matrix<double, 3, 3> mat;
using namespace Eigen;
using namespace std;

class state {
public:

	// state components
	Eigen::Matrix<double, 3, 3> M1, M2;

	// constructors
	state() : M1(), M2() {};      // constructors
	state(mat M1in, mat M2in) : M1(M1in), M2(M2in) {};


	// in place addition and multiplication
	state operator+=(const state & X) {
		M1 += X.M1; M2 += X.M2;
		return *this;
	}

	state operator*=(const double a) {
		M1 *= a; M2 *= a;
		return *this;
	}

	state operator+(const state &lhs, double rhs) {
		return state(lhs.M1 + rhs, lhs.M2 + rhs);
	}

	state operator+(double lhs, const state &rhs) {
		return state(lhs + rhs.M1, lhs + rhs.M2);
	}

	// ODE
	void operator() (const state & X, state & dX, const double) {
		dX.M1 = X.M1*X.M2.adjoint()*X.M2;
		dX.M2 = X.M2*X.M1.adjoint()*X.M1;
	}
};


// vector space operations

state operator+(const state &lhs, const state &rhs) {
	return state(lhs.M1 + rhs.M1, lhs.M2 + rhs.M2);
}

state operator*(const state &lhs, const double &rhs) {
	return state(lhs.M1*rhs, lhs.M2*rhs);
}

state operator*(const double &lhs, const state &rhs) {
	return state(lhs*rhs.M1, lhs*rhs.M2);
}

state operator/(const state &lhs, const state &rhs) {
	return state(lhs.M1.cwiseQuotient(rhs.M1), lhs.M2.cwiseQuotient(rhs.M2));
}
state abs(const state &X) {
	return state(abs(X.M1), abs(X.M2));
}


// lp infinity norm
namespace boost::numeric::odeint {
	template<>
	struct vector_space_norm_inf< state > {
		typedef double result_type;
		double operator()(const state &X) const {
			return max(X.M1.lpNorm<Infinity>(), X.M2.lpNorm<Infinity>());
		}
	};
}


//write to std output
void write(state &x, const double t) {
	cout << t << "\t" << x.M1 << "\t" << x.M2 << "\n";
}


//
// int main
//
int main(int argc, char* argv[]) {

	// set values
	mat M1, M2;
	double t_end = 1;
	double t_start = 10;

	M1 << 0.1, 0, 0, 0, 0.2, 0.1, 0.2, 0, 0.3;
	M2 << 0.5, 0, 0, 0, 0.6, 0, 0, 0, 0.7;

	state values(M1, M2);

	using namespace boost::numeric::odeint;

	// type definition for numerical integration
	typedef runge_kutta_dopri5< state, double, state, double, vector_space_algebra > stepper;

	// integration
	int steps = integrate_adaptive(make_controlled<stepper>(1E-10, 1E-10), state(), values, t_start, t_end, 0.01);

	//output
	write(values, t_end);
	return(0);
}