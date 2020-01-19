import pathlib
from common.file_system import get_root
from common.config import get_global_config,  get_odeint_config
import os.path
import numpy as np
import math

segment = 'medium'

system = 'dimer'
task = 'odeint_rk4'
is_continue = 'true'

step = (2.0 * math.pi) * 0.0005
total_num_periods = 1000

current_num_periods = 1000
current_num_time_points = 1001

Us = list(np.linspace(0.01, 1.0, 100, dtype=float))
Ns = list(np.linspace(100, 100, 1, dtype=int))

diss_type = 1
diss_gamma = 0.1

E = 0.0
J = 1.0

drv_type = 1
drv_ampl = 3.4
drv_freq = 1.0
drv_phase = 0.0

for N in Ns:
    for U in Us:

        print("N = " + str(N))
        print("U = " + str(U))
        local_path = '/' + system
        if task == 'odeint_rk4':
            local_path += '/' + task + '_' + str(total_num_periods) + '_' + str(format(step, '0.2e'))
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
        config_list += get_odeint_config(step, total_num_periods, current_num_periods, is_continue, data_path + '/')

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

