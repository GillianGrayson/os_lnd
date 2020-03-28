import pathlib
from common.file_system import get_root
from common.config import get_global_config,  get_odeint_config, get_smallest_eigen_vector_config, get_all_evals_config
import os.path
import numpy as np
import math

segment = 'medium'

system = 'dimer'
task = 'smallest_eigen_vector'
is_continue = 'true'

step = (2.0 * math.pi) * 0.0005
total_num_periods = 100

current_num_periods = 100
current_num_time_points = 101

max_num_iterations = 500000
tolerance = 1e-10

Us = list(np.linspace(1.0, 1.0, 1, dtype=float))
Ns = list(np.linspace(100, 100, 1, dtype=int))

diss_type = 1
diss_gamma = 0.1

E = 1.0
J = 1.0

drv_type = 0
drv_ampl = 1.5
drv_freq = 1.0
drv_phase = 0.0

for N in Ns:
    for U in Us:

        print("N = " + str(N))
        print("U = " + str(U))
        local_path = '/' + system
        if task == 'odeint':
            local_path += '/' + task + '_' + str(total_num_periods) + '_' + str(format(step, '0.2e'))
        elif task == 'smallest_eigen_vector':
            local_path += '/' + task + '_' + str(max_num_iterations) + '_' + str(format(tolerance, '0.2e'))
        elif task == 'all_evals':
            local_path += '/' + task + '_' + str(max_num_iterations) + '_' + str(format(tolerance, '0.2e'))
        else:
            local_path += '/' + task

        local_path += \
            '/np_' + str(N) + \
            '/diss_' + str(diss_type) + '_' + str(format(diss_gamma, '0.4f')) + \
            '/prm_' + str(format(E, '0.4f')) + '_' + str(format(U, '0.4f')) + '_' + str(format(J, '0.4f')) + \
            '/drv_' + str(drv_type) + '_' + str(format(drv_ampl, '0.4f')) + '_' + str(format(drv_freq, '0.4f')) + '_' + str(format(drv_phase, '0.4f'))

        data_path = get_root() + local_path

        config_list = []
        config_list.append('[dimer]')
        config_list.append('num_particles = ' + str(N))
        config_list.append('diss_type = ' + str(diss_type))
        config_list.append('diss_gamma = ' + str(diss_gamma))
        config_list.append('E = ' + str(E))
        config_list.append('U = ' + str(U))
        config_list.append('J = ' + str(J))
        config_list.append('drv_type = ' + str(drv_type))
        config_list.append('drv_ampl = ' + str(drv_ampl))
        config_list.append('drv_freq = ' + str(drv_freq))
        config_list.append('drv_phase = ' + str(drv_phase))

        config_list += get_global_config(system, task)
        if task == 'odeint':
            config_list += get_odeint_config(step, total_num_periods, current_num_periods, current_num_time_points, is_continue, data_path + '/')
        elif task == 'smallest_eigen_vector':
            config_list += get_smallest_eigen_vector_config(max_num_iterations, tolerance)
        elif task == 'all_evals':
            config_list += get_all_evals_config(max_num_iterations, tolerance)

        pathlib.Path(data_path).mkdir(parents=True, exist_ok=True)

        file_config = open(data_path + '/config.ini', 'w')
        file_config.write('\n'.join(config_list))

        fn_suffix = \
            'np(' + str(N) + ')_' + \
            'diss(' + str(diss_type) + '_' + str(format(diss_gamma, '0.4f')) + ')_' + \
            'prm(' + str(format(E, '0.4f')) + '_' + str(format(U, '0.4f')) + '_' + str(format(J, '0.4f')) + ')_' + \
            'drv(' + str(drv_type) + '_' + str(format(drv_ampl, '0.4f')) + '_' + str(
                format(drv_freq, '0.4f')) + '_' + str(format(drv_phase, '0.4f')) + ')'

        fn_test = data_path + '/rho_mtx_' + fn_suffix + '.txt'

        if not os.path.isfile(fn_test) or is_continue == 'true':
            if segment == 'short':
                os.system('sbatch run_short.sh ' + data_path)
            elif segment == 'medium':
                os.system('sbatch run_medium.sh ' + data_path)

