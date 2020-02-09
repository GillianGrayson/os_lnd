#define FMT_HEADER_ONLY
#define FMT_DEPRECATED
#define EIGEN_USE_MKL_ALL

#include "setup_processor.h"
#include "run_processor.h"

#include <slepceps.h>

int main(int argc, char* argv[])
{
	PetscErrorCode ierr;
	static char help[] = "";
	ierr = SlepcInitialize(&argc, &argv, (char*)0, help); if (ierr) return ierr;
	
	INIReader ini("config.ini");
	if (ini.ParseError() < 0)
	{
		throw std::runtime_error("Can't load 'test.ini'\n");
	}

	Model model(ini);

	SetupProcessor setup_processor;
	setup_processor.set_strategy(model);
	setup_processor.process(model);

	RunProcessor run_processor;
	run_processor.set_strategy(model);
	run_processor.process(model);

	ierr = SlepcFinalize();
}