import pathlib
from common.file_system import get_root
from common.config import get_regular_global_config
import os.path
import numpy as np
import math

segment = 'medium'

system = 'super_decoh'

task = 'eigen_dense'

method = 'simple'
G_type = 0
aux_dim = 0
reshuffle_type = 0

#ps = list(np.logspace(-3.0, 0.0, num=31, base=10.0))
ps = [0.8]
#Ns = list(np.logspace(1, 2, 11, base=10.0, dtype=int))
Ns = [100]
total_num_evals  = 100000
#all_seeds = [list(np.linspace(1, math.ceil(total_num_evals/(N * N)), math.ceil(total_num_evals/(N * N)), dtype=int))  for N in Ns]
all_seeds = [list(np.linspace(1, 1000, 1000, dtype=int))]

#ps = [1.0]
#Ns = [200]
#all_seeds = [list(np.linspace(1, 20, 20, dtype=int))]

num_seeds = 1000000

for N_id, N in enumerate(Ns):
    seeds = all_seeds[N_id]
    for p in ps:
        for seed in seeds:

            print("N = " + str(N))
            print("p = " + str(p))
            print("seed = " + str(seed))
            local_path = '/regular/' + system
            local_path += '/' + task

            local_path += \
                '/method_' + str(method) + \
                '/G_type_' + str(G_type) + '_ad_' + str(aux_dim) + \
                '/reshuffle_type_' + str(reshuffle_type) + \
                '/N_' + str(N) + \
                '/p_' + str(format(p, '0.10f')) + \
                '/seed_' + str(seed)

            data_path = get_root() + local_path

            config_list = []
            config_list.append('[super_decoh]')
            config_list.append('method = ' + str(method))
            config_list.append('N = ' + str(N))
            config_list.append('aux_dim = ' + str(aux_dim))
            config_list.append('seed = ' + str(seed))
            config_list.append('num_seeds = ' + str(num_seeds))
            config_list.append('p = ' + str(p))
            config_list.append('G_type = ' + str(G_type))
            config_list.append('reshuffle_type = ' + str(reshuffle_type))
            config_list.append('save_G = false')
            config_list.append('evals_G = false')
            config_list.append('save_A = false')

            config_list += get_regular_global_config(system, task, save_rho='false', name_precision=10)

            pathlib.Path(data_path).mkdir(parents=True, exist_ok=True)

            file_config = open(data_path + '/config.ini', 'w')
            file_config.write('\n'.join(config_list))

            fn_suffix = \
                'reshuffle(' + str(reshuffle_type) + ')_' + \
                'G(' + str(G_type) + ')_' + \
                'N(' + str(N) + ')_' + \
                'ad(' + str(aux_dim) + ')_' + \
                'p('  + str(format(p, '0.10f')) + ')_' + \
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
