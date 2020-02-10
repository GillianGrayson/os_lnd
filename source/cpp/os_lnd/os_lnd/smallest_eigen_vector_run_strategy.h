#pragma once
#include "run_strategy.h"
#include "save.h"

#include <slepceps.h>

struct SmallestEigenVectorRunStrategy : RunStrategy
{
	void run(Model& model) override
	{
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		
        Mat            A;           /* problem matrix */
		MatInfo		   mat_info;
        EPS            eps;         /* eigenproblem solver context */
        EPSType        type;
        PetscReal      error, tol, re, im;
        PetscScalar    value, kr, ki;
		PetscScalar    *sm_vector;
        Vec            xr, xi;
        PetscInt       n, i, Istart, Iend, nev, ncv, maxit, its, nconv, size_n, size_m;

		n = model.sys_size * model.sys_size;

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		   Compute the operator matrix that defines the eigensystem, Ax=kx
		   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        MatCreate(PETSC_COMM_WORLD, &A);
        MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
        MatSetFromOptions(A);
        MatSetUp(A);
		for (int k = 0; k < model.lindbladian.outerSize(); ++k)
		{
			for (typename sp_mtx::InnerIterator it(model.lindbladian, k); it; ++it)
			{
				value = it.value().real() + it.value().imag() * PETSC_i;			
				MatSetValue(A, it.row(), it.col(), value, INSERT_VALUES);
				//model.log_message(fmt::format("a[{:d},{:d}] = {:16e} + {:16e} i", it.row(), it.col(), PetscRealPart(value), PetscImaginaryPart(value)));
			}
		}
		MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

		MatView(A, PETSC_VIEWER_STDOUT_WORLD);

		MatGetInfo(A, MAT_GLOBAL_MAX, &mat_info);
		model.log_message(fmt::format("mallocs: {:16e}", mat_info.mallocs));
		model.log_message(fmt::format("nz_allocated: {:16e}", mat_info.nz_allocated));
		model.log_message(fmt::format("nz_used: {:16e}", mat_info.nz_used));

		MatGetSize(A, &size_n, &size_m);
		model.log_message(fmt::format("size_n: {:d}", size_n));
		model.log_message(fmt::format("size_m: {:d}", size_m));

		MatCreateVecs(A, NULL, &xr);
		MatCreateVecs(A, NULL, &xi);

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
					  Create the eigensolver and set various options
		   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
		EPSCreate(PETSC_COMM_WORLD, &eps);
		EPSSetOperators(eps, A, NULL);
		EPSSetProblemType(eps, EPS_NHEP);
		EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);
		EPSSetType(eps, EPSKRYLOVSCHUR);
		EPSSetTolerances(eps, 1e-8, 500000);
		EPSSetFromOptions(eps);


		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
							Solve the eigensystem
		   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
		EPSSolve(eps);
		/*
		   Optional: Get some information from the solver and display it
		*/
		EPSGetIterationNumber(eps, &its);
		model.log_message(fmt::format("Number of iterations of the method: {:d}", its));
		EPSGetType(eps, &type);
		model.log_message(fmt::format("Solution method: {:s}", type));
		EPSGetDimensions(eps, &nev, &ncv, NULL);
		model.log_message(fmt::format("Number of requested eigenvalues: {:d}", nev));
		model.log_message(fmt::format("Maximum dimension of the subspace to be used by the solver: {:d}", ncv));
		EPSGetTolerances(eps, &tol, &maxit);
		model.log_message(fmt::format("Stopping condition: tol{:16e}", (double)tol));
		model.log_message(fmt::format("Stopping condition: maxit={:d}", maxit));

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
							Display solution and clean up
		   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
		/*
				Get number of converged approximate eigenpairs
		*/
		EPSGetConverged(eps, &nconv);
		model.log_message(fmt::format("Number of converged eigenpairs: {:d}", nconv));

		/*
		   Display eigenvalue and relative error
		*/

		EPSGetEigenpair(eps, 0, &kr, &ki, xr, xi);
		/*
		   Compute the relative error associated to eigenpair
		*/
		EPSComputeError(eps, 0, EPS_ERROR_RELATIVE, &error);

		re = PetscRealPart(kr);
		im = PetscImaginaryPart(kr);

		if (im != 0.0) {
			PetscPrintf(PETSC_COMM_WORLD, " %9f%+9fi %12g\n", (double)re, (double)im, (double)error);
		}
		else {
			PetscPrintf(PETSC_COMM_WORLD, "   %12f       %12g\n", (double)re, (double)error);
		}
		
		VecGetArray(xr, &sm_vector);

		Eigen::VectorXcd rho_vec(model.sys_size * model.sys_size);
		for (int st_id = 0; st_id < model.sys_size * model.sys_size; st_id++)
		{
			double real_part = double(PetscRealPart(sm_vector[st_id]));
			double imag_part = double(PetscImaginaryPart(sm_vector[st_id]));
			rho_vec[st_id] = std::complex<double>(real_part, imag_part);
		}
		
		VecRestoreArray(xr, &sm_vector);

		model.rho = Eigen::Map<Eigen::MatrixXcd>(rho_vec.data(), model.sys_size, model.sys_size);
		auto trace_rho = model.rho.trace();
		model.rho = model.rho / trace_rho;

		auto fn = "rho_mtx" + model.suffix;
		save_dense_mtx(model.rho, fn, save_precision);
		
		model.log_time_duration();
		model.log_memory_usage();

		/*
		   Free work space
		*/
		EPSDestroy(&eps);
		MatDestroy(&A);
		VecDestroy(&xr);
		VecDestroy(&xi);
	}
};