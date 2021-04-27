import pathlib
from common.file_system import get_root
from common.config import get_serial_global_config, get_odeint_config
import os.path
import numpy as np
import math

segment = 'short'

system = 'xxz'

task = 'odeint'

silent = 'false'

num_spins = 7
mu = 0.001
drv_type = 0
ampl = 1.0
freq = 2 * np.pi
phase = 0.0
Deltas = list(np.linspace(0.0, 2.0, 51, dtype=float))
Ws = list(np.linspace(0.0, 1.0, 51, dtype=float))

T = 2.0 * np.pi / freq

start_state_id = 0
step = T / 25
dump_type = 'linear'
num_trans_periods = 500
num_obser_periods = 1
current_num_obser_periods = 1
current_num_obser_time_points = 21
dump_progress = 'false'
dump_last_time = 'false'
is_continue = 'false'
continue_path = ''

serial_start = 1
serial_shift = 1
serial_num = 100
num_seeds = 1000000


for W in Ws:
    for Delta in Deltas:
        print("W = " + str(W))
        print("Delta = " + str(Delta))

        local_path = '/serial/' + system
        local_path += '/' + task

        local_path += \
            '/n_' + str(num_spins)  + \
            '/params_' + str(format(Delta, '0.4f')) + '_' + str(format(W, '0.4f')) + '_' + str(format(mu, '0.4f')) + '_' + str(drv_type) + '_' + str(format(ampl, '0.4f')) + '_' + str(format(freq, '0.4f')) + '_' + str(format(phase, '0.4f')) + \
            '/seeds_' + str(serial_start) + '_' + str(serial_shift) + '_' + str(serial_num)

        data_path = get_root() + local_path

        config_list = []
        config_list.append('[xxz]')
        config_list.append('num_spins = ' + str(num_spins))
        config_list.append('seed = ' + str(serial_start))
        config_list.append('num_seeds = ' + str(num_seeds))
        config_list.append('Delta = ' + str(Delta))
        config_list.append('W = ' + str(W))
        config_list.append('mu = ' + str(mu))
        config_list.append('drv_type = ' + str(drv_type))
        config_list.append('ampl = ' + str(ampl))
        config_list.append('freq = ' + str(freq))
        config_list.append('phase = ' + str(phase))

        config_list += get_serial_global_config(system, task, serial_start, serial_shift, serial_num, name_precision=4, silent=silent)

        config_list += get_odeint_config(
            step,
            num_obser_periods,
            num_trans_periods,
            current_num_obser_periods,
            current_num_obser_time_points,
            is_continue,
            continue_path,
            dump_progress,
            dump_last_time
        )

        pathlib.Path(data_path).mkdir(parents=True, exist_ok=True)

        file_config = open(data_path + '/config.ini', 'w')
        file_config.write('\n'.join(config_list))

        fn_suffix = \
            'serial(' + str(format(serial_start, '0.4f')) + '_' + str(format(serial_shift, '0.4f')) + '_' + str(serial_num) + ')_' + \
            'ns(' + str(num_spins) + ')_' + \
            'prm(' + str(format(Delta, '0.4f')) + '_' + str(format(W, '0.4f')) + '_' + str(format(mu, '0.4f')) + '_' + str(drv_type) + '_' + str(format(ampl, '0.4f')) + '_' + str(format(freq, '0.4f')) + '_' + str(format(phase, '0.4f')) + ')'

        fn_test = data_path + '/serial_j_0_' + fn_suffix + '.txt'

        if not os.path.isfile(fn_test):
            # print(fn_test)
            if segment == 'short':
                os.system('sbatch run_short.sh ' + data_path)
            elif segment == 'medium':
                os.system('sbatch run_medium.sh ' + data_path)
            elif segment == 'long':
                os.system('sbatch run_long.sh ' + data_path)
