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
        EPS            eps;         /* eigenproblem solver context */
        EPSType        type;
        PetscReal      error, tol, re, im;
        PetscScalar    value, kr, ki;
        Vec            xr, xi;
        PetscInt       n, i, Istart, Iend, nev, maxit, its, nconv;
        PetscErrorCode ierr;

		n = model.sys_size;

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		   Compute the operator matrix that defines the eigensystem, Ax=kx
		   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
        ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n); CHKERRQ(ierr);
        ierr = MatSetFromOptions(A); CHKERRQ(ierr);
        ierr = MatSetUp(A); CHKERRQ(ierr);
		for (int k = 0; k < model.lindbladian.outerSize(); ++k)
		{
			for (typename sp_mtx::InnerIterator it(model.lindbladian, k); it; ++it)
			{
				value = it.value().real() + it.value().imag() * PETSC_i;			
				ierr = MatSetValue(A, it.row(), it.col(), value, INSERT_VALUES); CHKERRQ(ierr);
			}
		}
		ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
		ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

		ierr = MatCreateVecs(A, NULL, &xr); CHKERRQ(ierr);
		ierr = MatCreateVecs(A, NULL, &xi); CHKERRQ(ierr);

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
					  Create the eigensolver and set various options
		   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
		ierr = EPSCreate(PETSC_COMM_WORLD, &eps); CHKERRQ(ierr);
		ierr = EPSSetOperators(eps, A, NULL); CHKERRQ(ierr);
		ierr = EPSSetProblemType(eps, EPS_NHEP); CHKERRQ(ierr);
		ierr = EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE); CHKERRQ(ierr);
		ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);


		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
							Solve the eigensystem
		   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
		ierr = EPSSolve(eps); CHKERRQ(ierr);
		/*
		   Optional: Get some information from the solver and display it
		*/
		ierr = EPSGetIterationNumber(eps, &its); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of iterations of the method: %D\n", its); CHKERRQ(ierr);
		ierr = EPSGetType(eps, &type); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD, " Solution method: %s\n\n", type); CHKERRQ(ierr);
		ierr = EPSGetDimensions(eps, &nev, NULL, NULL); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of requested eigenvalues: %D\n", nev); CHKERRQ(ierr);
		ierr = EPSGetTolerances(eps, &tol, &maxit); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD, " Stopping condition: tol=%.4g, maxit=%D\n", (double)tol, maxit); CHKERRQ(ierr);

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
							Display solution and clean up
		   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
		/*
				Get number of converged approximate eigenpairs
		*/
		ierr = EPSGetConverged(eps, &nconv); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of converged eigenpairs: %D\n\n", nconv); CHKERRQ(ierr);

		if (nconv > 0) 
		{
			/*
			   Display eigenvalue and relative error
			*/

			ierr = EPSGetEigenpair(eps, 0, &kr, &ki, xr, xi); CHKERRQ(ierr);
			/*
			   Compute the relative error associated to eigenpair
			*/
			ierr = EPSComputeError(eps, 0, EPS_ERROR_RELATIVE, &error); CHKERRQ(ierr);

			re = PetscRealPart(kr);
			im = PetscImaginaryPart(kr);

			if (im != 0.0) {
				ierr = PetscPrintf(PETSC_COMM_WORLD, " %9f%+9fi %12g\n", (double)re, (double)im, (double)error); CHKERRQ(ierr);
			}
			else {
				ierr = PetscPrintf(PETSC_COMM_WORLD, "   %12f       %12g\n", (double)re, (double)error); CHKERRQ(ierr);
			}
		}

		Eigen::VectorXcd rho_vec(model.sys_size * model.sys_size);
		for (int st_id = 0; st_id < model.sys_size * model.sys_size; st_id++)
		{
			rho_vec[st_id] = std::complex<double>(PetscRealPart(xr[st_id]), PetscImaginaryPart(xr[st_id]));
		}

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
		ierr = EPSDestroy(&eps); CHKERRQ(ierr);
		ierr = MatDestroy(&A); CHKERRQ(ierr);
		ierr = VecDestroy(&xr); CHKERRQ(ierr);
		ierr = VecDestroy(&xi); CHKERRQ(ierr);
	}
};