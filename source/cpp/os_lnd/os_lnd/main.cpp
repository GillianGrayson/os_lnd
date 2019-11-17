#define FMT_HEADER_ONLY
#define FMT_DEPRECATED
#define ARMA_DONT_USE_OPENMP 1
#define ARMA_USE_ARPACK 1

#include "setup_processor.h"
#include "run_processor.h"

int main(int argc, char* argv[])
{
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
}