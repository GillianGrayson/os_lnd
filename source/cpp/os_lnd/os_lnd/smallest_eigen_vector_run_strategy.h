#pragma once
#include "run_strategy.h"
#include "save.h"

#ifdef __linux__ 
#include <slepceps.h>
#endif

struct SmallestEigenVectorRunStrategy : RunStrategy
{
	void run(Model& model) override
	{
#ifdef __linux__
		const int save_precision = model.ini.GetInteger("global", "save_precision", 0);
		const int max_num_iterations = model.ini.GetInteger("smallest_eigen_vector", "max_num_iterations", 0);
		const double tolerance = model.ini.GetReal("smallest_eigen_vector", "tolerance", 0.0);
		
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
		model.log_message(fmt::format("SLEPc matrix init ..."));
		model.log_time_duration();
        MatCreate(PETSC_COMM_WORLD, &A);
        MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
        MatSetFromOptions(A);
        MatSetUp(A);
		for (int k = 0; k < model.lindbladian.outerSize(); ++k)
		{
			std::vector<std::complex<double>> values;
			std::vector<int> rows;
			for (typename sp_mtx::InnerIterator it(model.lindbladian, k); it; ++it)
			{
				values.push_back(it.value());
				rows.push_back(it.row());
				
				//model.log_message(fmt::format("a[{:d},{:d}] = {:16e} + {:16e} i", it.row(), it.col(), it.value().real(), it.value().imag()));
				//model.log_message(fmt::format("a[{:d},{:d}] = {:16e} + {:16e} i", it.row(), it.col(), PetscRealPart(value), PetscImaginaryPart(value)));
			}

			int num_rows = rows.size();

			PetscScalar* values_petsc = new PetscScalar[num_cols];
			PetscInt* rows_petsc = new PetscInt[num_cols];
			for (int r_id = 0; r_id < num_rows; r_id++)
			{
				rows_petsc[r_id] = rows[r_id];
				values_petsc[r_id] = values[r_id].real() + values[r_id].imag() * PETSC_i;

				//model.log_message(fmt::format("ba[{:d}] = {:16e} + {:16e} i", cols[c_id], PetscRealPart(values_petsc[c_id]), PetscImaginaryPart(values_petsc[c_id])));
			}
			MatSetValues(A, 1, &k, num_rows, rows_petsc, values_petsc, INSERT_VALUES);

			delete[] values_petsc;
			delete[] rows_petsc;

			//for (typename sp_mtx::InnerIterator it(model.lindbladian, k); it; ++it)
			//{
			//	value = it.value().real() + it.value().imag() * PETSC_i;			
			//	MatSetValue(A, it.row(), it.col(), value, INSERT_VALUES);
			//	//model.log_message(fmt::format("a[{:d},{:d}] = {:16e} + {:16e} i", it.row(), it.col(), PetscRealPart(value), PetscImaginaryPart(value)));
			//}
		}
		MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
		model.log_message(fmt::format("SLEPc matrix init end"));
		model.log_time_duration();

		//MatView(A, PETSC_VIEWER_STDOUT_WORLD);

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
		EPSSetTolerances(eps, tolerance, max_num_iterations);
		EPSSetFromOptions(eps);
		
		EPSGetTolerances(eps, &tol, &maxit);
		model.log_message(fmt::format("Stopping condition: tolerance = {:16e}", (double)tol));
		model.log_message(fmt::format("Stopping condition: max_num_iterations = {:d}", maxit));

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

		model.log_message(fmt::format("smallest_eval = {:16e} + {:16e} i", re, im));
		model.log_message(fmt::format("error = {:16e}", error));
		
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
		std::complex<double> trace_rho = model.rho.trace();
		model.log_message(fmt::format("trace_rho = {:16e} + {:16e} i", trace_rho.real(), trace_rho.imag()));
		model.rho = model.rho / trace_rho;
		trace_rho = model.rho.trace();
		model.log_message(fmt::format("trace_rho = {:16e} + {:16e} i", trace_rho.real(), trace_rho.imag()));

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

		ModelProcessor model_processor;
		model_processor.set_strategy(model);
		model_processor.init_model(model);
		model_processor.release_observables(model);
		model.log_time_duration();
#endif
	}
};