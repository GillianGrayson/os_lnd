#include "model_setup_processor.h"
#include "model_run_processor.h"


int main(int argc, char* argv[])
{
	INIReader ini("config.ini");
	if (ini.ParseError() < 0)
	{
		throw std::runtime_error("Can't load 'test.ini'\n");
	}

	Model model(ini);

	ModelSetupProcessor setup_processor;
	setup_processor.set_model_strategy(model);
	setup_processor.process_model(model);

	ModelRunProcessor run_processor;
	run_processor.set_model_strategy(model);
	run_processor.process_model(model);
}