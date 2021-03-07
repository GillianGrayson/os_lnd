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

	std::string logger_type = ini.Get("global", "logger_type", "unknown");
	auto logger = spdlog::stdout_color_mt(logger_type);

	const int save_precision = ini.GetInteger("global", "save_precision", 0);
	const int name_precision = ini.GetInteger("global", "name_precision", 0);

	const std::string run_type = ini.Get("global", "run_type", "unknown");

	const auto serial_start = ini.GetReal("global", "serial_start", 0.0);
	const auto serial_shift = ini.GetReal("global", "serial_shift", 0.0);
	const int serial_num = ini.GetInteger("global", "serial_num", 1);

	if (run_type == "regular")
	{
		Model model(ini);

		ModelProcessor setup_processor;
		setup_processor.set_strategy(model);
		setup_processor.create_model(model);

		RunProcessor run_processor;
		run_processor.set_strategy(model);
		run_processor.process(model);
	}
	else if (run_type == "serial")
	{
		std::string format = fmt::format("0:.{:d}f", name_precision);
		auto suffix = fmt::format(fmt::format("_serial({{{:s}}}", format), serial_start) + fmt::format(fmt::format("_{{{:s}}}", format), serial_shift) + fmt::format("_{:d})", serial_num);

		std::map<std::string, std::vector<double>> features_double;
		std::map<std::string, std::vector<std::complex<double>>> features_complex;

		for (int serial_id = 0; serial_id < serial_num; serial_id++)
		{
			Model model(ini);

			ModelProcessor setup_processor;
			setup_processor.set_strategy(model);

			RunProcessor run_processor;
			run_processor.set_strategy(model);
			
			double serial_state = serial_start + serial_shift * serial_id;
			model.set_serial_state(serial_state);
			setup_processor.create_model(model);
			
			if (serial_id == 0)
			{
				suffix += model.suffix_serial;
				setup_processor.set_serial_features(model, features_double, features_complex);
			}
			
			run_processor.process_serial(model, features_double, features_complex);
		}

		for (auto const& x : features_double)
		{
			auto fn = "serial_" + x.first + suffix;
			save_vector(x.second, fn, save_precision);
		}

		for (auto const& x : features_complex)
		{
			auto fn = "serial_" + x.first + suffix;
			save_vector(x.second, fn, save_precision);
		}
	}
	else
	{
		throw std::runtime_error("Unsupported run_type");
	}

	
#ifdef __linux__ 
	ierr = SlepcFinalize();
#endif
}