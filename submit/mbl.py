import pathlib
from common.file_system import get_root
from common.config import get_global_config,  get_odeint_config, get_smallest_eigen_vector_config
import os.path
import numpy as np
import math

segment = 'medium'

system = 'mbl'

task = 'smallest_eigen_vector'

is_continue = 'true'
step = 0.01
total_num_periods = 1000

current_num_periods = 1000

max_num_iterations = 500000
tolerance = 1e-10

num_spins = 8

seeds = list(np.linspace(10, 10, 1, dtype=int))
num_seeds = 1000000

diss_type = 1
diss_phase = 0.0
diss_gamma = 0.1

Ws = list(np.linspace(2.0, 2.0, 1, dtype=float))
U = 1.0
J = 1.0

for seed in seeds:
    for W in Ws:

        print("seed = " + str(seed))
        print("W = " + str(W))
        local_path = '/' + system
        if task == 'odeint_rk4':
            local_path += '/' + task + '_' + str(total_num_periods) + '_' + str(format(step, '0.2e'))
        elif task == 'smallest_eigen_vector':
            local_path += '/' + task + '_' + str(max_num_iterations) + '_' + str(format(tolerance, '0.2e'))
        else:
            local_path += '/' + task

        local_path += \
            '/ns_' + str(num_spins) + \
            '/diss_' + str(diss_type) + '_' + str(format(diss_phase, '0.4f')) + '_' + str(format(diss_gamma, '0.4f')) + \
            '/prm_' + str(format(W, '0.4f')) + '_' + str(format(U, '0.4f')) + '_' + str(format(J, '0.4f')) + \
            '/seed_' + str(seed)

        data_path = get_root() + local_path

        config_list = []
        config_list.append('[mbl]')
        config_list.append('num_spins = ' + str(num_spins))
        config_list.append('seed = ' + str(seed))
        config_list.append('num_seeds = ' + str(num_seeds))
        config_list.append('diss_type = ' + str(diss_type))
        config_list.append('diss_phase = ' + str(diss_phase))
        config_list.append('diss_gamma = ' + str(diss_gamma))
        config_list.append('W = ' + str(W))
        config_list.append('U = ' + str(U))
        config_list.append('J = ' + str(J))

        config_list += get_global_config(system, task)
        if task == 'odeint_rk4':
            config_list += get_odeint_config(step, total_num_periods, current_num_periods, is_continue, data_path + '/')
        elif task == 'smallest_eigen_vector':
            config_list += get_smallest_eigen_vector_config(max_num_iterations, tolerance)

        pathlib.Path(data_path).mkdir(parents=True, exist_ok=True)

        file_config = open(data_path + '/config.ini', 'w')
        file_config.write('\n'.join(config_list))

        fn_suffix = \
            'ns(' + str(num_spins) + ')_' + \
            'seed(' + str(seed) + ')_' + \
            'diss(' + str(diss_type) + '_' + str(format(diss_phase, '0.4f')) + '_' + str(format(diss_gamma, '0.4f')) + ')_' + \
            'prm(' + str(format(W, '0.4f')) + '_' + str(format(U, '0.4f')) + '_' + str(format(J, '0.4f')) + ')'

        fn_test = data_path + '/rho_mtx_' + fn_suffix + '.txt'

        if not os.path.isfile(fn_test) or is_continue == 'true':
            if segment == 'short':
                os.system('sbatch run_short.sh ' + data_path)
            elif segment == 'medium':
                os.system('sbatch run_medium.sh ' + data_path)

