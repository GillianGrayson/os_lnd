#include "model_processor.h"


int main(int argc, char* argv[])
{
	INIReader ini("config.ini");
	if (ini.ParseError() < 0)
	{
		throw std::runtime_error("Can't load 'test.ini'\n");
	}

	Model model(ini);

	ModelProcessor model_processor;
	model_processor.set_model_strategy(model);
	model_processor.process_model(model);
}
