import pathlib
from common.file_system import get_root
from common.config import get_regular_global_config
import os.path
import numpy as np
import math

segment = 'short'

system = 'lfk'

task = 'eigen_dense'

Ns = [100]
total_num_evals  = 100000
all_seeds = [list(np.linspace(1, 1000, 1000, dtype=int))]

num_seeds = 1000000

for N_id, N in enumerate(Ns):
    seeds = all_seeds[N_id]
    for seed in seeds:

        print("N = " + str(N))
        print("seed = " + str(seed))
        local_path = '/regular/' + system
        local_path += '/' + task

        local_path += \
            '/N_' + str(N) + \
            '/seed_' + str(seed)

        data_path = get_root() + local_path

        config_list = []
        config_list.append('[lfk]')
        config_list.append('N = ' + str(N))
        config_list.append('seed = ' + str(seed))
        config_list.append('num_seeds = ' + str(num_seeds))

        config_list += get_regular_global_config(system, task, save_rho='false', name_precision=10)

        pathlib.Path(data_path).mkdir(parents=True, exist_ok=True)

        file_config = open(data_path + '/config.ini', 'w')
        file_config.write('\n'.join(config_list))

        fn_suffix = \
            'N(' + str(N) + ')_' + \
            'seed(' + str(seed) + ')'

        fn_test = data_path + '/lindbladian_evals_' + fn_suffix + '.txt'

        if not os.path.isfile(fn_test):
            print('file does not exist')
            if segment == 'short':
                os.system('sbatch run_short.sh ' + data_path)
            elif segment == 'medium':
                os.system('sbatch run_medium.sh ' + data_path)
            elif segment == 'long':
                os.system('sbatch run_long.sh ' + data_path)
