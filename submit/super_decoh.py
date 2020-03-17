import pathlib
from common.file_system import get_root
from common.config import get_global_config
import os.path
import numpy as np

segment = 'medium'

system = 'super_decoh'

task = 'eigen_dense'

reshuffle_type = 0
Ns = list(np.linspace(10, 100, 91, dtype=int))
#ps = list(np.linspace(0.1, 1.0, 10, dtype=float))
#ps = list(np.logspace(-10.0, 0.0, num=11, base=10.0))
ps = [1e-10]
seeds = list(np.linspace(1, 100, 100, dtype=int))
num_seeds = 1000000

for N in Ns:
    for p in ps:
        for seed in seeds:

            print("N = " + str(N))
            print("p = " + str(p))
            print("seed = " + str(seed))
            local_path = '/' + system
            local_path += '/' + task

            local_path += \
                '/reshuffle_type_' + str(reshuffle_type) + \
                '/N_' + str(N) + \
                '/p_' + str(format(p, '0.10f')) + \
                '/seed_' + str(seed)

            data_path = get_root() + local_path

            config_list = []
            config_list.append('[super_decoh]')
            config_list.append('N = ' + str(N))
            config_list.append('seed = ' + str(seed))
            config_list.append('num_seeds = ' + str(num_seeds))
            config_list.append('p = ' + str(p))
            config_list.append('reshuffle_type = ' + str(reshuffle_type))
            config_list.append('save_G = false')
            config_list.append('save_A = false')

            config_list += get_global_config(system, task, 10)

            pathlib.Path(data_path).mkdir(parents=True, exist_ok=True)

            file_config = open(data_path + '/config.ini', 'w')
            file_config.write('\n'.join(config_list))

            fn_suffix = \
                'reshuffle(' + str(reshuffle_type) + ')_' + \
                'N(' + str(N) + ')_' + \
                'p('  + str(format(p, '0.10f')) + ')_' + \
                'seed(' + str(seed) + ')'

            fn_test = data_path + '/rho_mtx_' + fn_suffix + '.txt'

            if not os.path.isfile(fn_test):
                if segment == 'short':
                    os.system('sbatch run_short.sh ' + data_path)
                elif segment == 'medium':
                    os.system('sbatch run_medium.sh ' + data_path)
