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

	const int save_precision = ini.GetInteger("global", "save_precision", 0);
	const int name_precision = ini.GetInteger("global", "name_precision", 0);

	const std::string run_type = ini.Get("global", "run_type", "unknown");

	const auto serial_start = ini.GetReal("global", "serial_start", 0.0);
	const auto serial_shift = ini.GetReal("global", "serial_shift", 0.0);
	const int serial_num = ini.GetInteger("global", "serial_num", 1);
	const auto serial_rho_evals = ini.GetBoolean("global", "serial_rho_evals", false);
	const auto serial_lindbladian_evals = ini.GetBoolean("global", "serial_lindbladian_evals", false);

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
		Model model(ini);
		
		std::string format = fmt::format("0:.{:d}f", name_precision);
		auto suffix = fmt::format(fmt::format("_serial({{{:s}}}", format), serial_start) + fmt::format(fmt::format("_{{{:s}}}", format), serial_shift) + fmt::format("_{:d})", serial_num);
		std::vector<std::complex<double>> rho_evals;
		std::vector<std::complex<double>> lindbladian_evals;
		
		for (int serial_id = 0; serial_id < serial_num; serial_id++)
		{
			double serial_state = serial_start + serial_shift * serial_id;

			model.set_serial_state(serial_state);

			ModelProcessor setup_processor;
			setup_processor.set_strategy(model);
			setup_processor.create_model(model);

			RunProcessor run_processor;
			run_processor.set_strategy(model);
			run_processor.process(model);

			if (serial_id == 0)
			{
				suffix += model.suffix_serial;
				
				if (serial_rho_evals)
				{
					rho_evals = model.rho_evals;
				}

				if (serial_lindbladian_evals)
				{
					lindbladian_evals = model.lindbladian_evals;
				}
			}
			else
			{
				if (serial_rho_evals)
				{
					rho_evals.insert(std::end(rho_evals), std::begin(model.rho_evals), std::end(model.rho_evals));
				}
				
				if (serial_lindbladian_evals)
				{
					lindbladian_evals.insert(std::end(lindbladian_evals), std::begin(model.lindbladian_evals), std::end(model.lindbladian_evals));
				}
			}
		}

		if (serial_rho_evals)
		{
			auto fn = "serial_rho_evals" + suffix;
			save_vector(rho_evals, fn, save_precision);
		}

		if (serial_lindbladian_evals)
		{
			auto fn = "serial_lindbladian_evals" + suffix;
			save_vector(lindbladian_evals, fn, save_precision);
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