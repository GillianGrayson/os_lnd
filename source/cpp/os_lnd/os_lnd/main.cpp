#define FMT_HEADER_ONLY
#define FMT_DEPRECATED
#define EIGEN_USE_MKL_ALL

#include "model_processor.h"
#include "run_processor.h"

#ifdef __linux__ 
#include <slepceps.h>
#endif

int main(int argc, char* argv[])
{
#ifdef __linux__ 
	PetscErrorCode ierr;
	static char help[] = "";
	ierr = SlepcInitialize(&argc, &argv, (char*)0, help); if (ierr) return ierr;
#endif
	
	INIReader ini("config.ini");
	if (ini.ParseError() < 0)
	{
		throw std::runtime_error("Can't load 'test.ini'\n");
	}

	Model model(ini);

	ModelProcessor setup_processor;
	setup_processor.set_strategy(model);
	setup_processor.create_model(model);

	RunProcessor run_processor;
	run_processor.set_strategy(model);
	run_processor.process(model);
	
#ifdef __linux__ 
	ierr = SlepcFinalize();
#endif
}