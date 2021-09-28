import pathlib
from common.file_system import get_root
from common.config import get_regular_global_config
import os.path
import numpy as np
import math

segment = 'medium'

system = 'lind_ham'

task = 'eigen_dense'

N = 100
alpha = 0.5
all_seeds = [list(np.linspace(1, 100, 100, dtype=int))]
num_seeds = 1000000

for seed in all_seeds:
    print("seed = " + str(seed))
    local_path = f"/regular/{system}/{task}/N_{N}/alpha_{alpha}/seed_{seed}"

    data_path = get_root() + local_path

    config_list = []
    config_list.append('[lind_ham]')
    config_list.append('N = ' + str(N))
    config_list.append('alpha = ' + str(alpha))
    config_list.append('seed = ' + str(seed))
    config_list.append('num_seeds = ' + str(num_seeds))

    config_list += get_regular_global_config(system, task, save_rho='true', name_precision=4)

    pathlib.Path(data_path).mkdir(parents=True, exist_ok=True)

    file_config = open(data_path + '/config.ini', 'w')
    file_config.write('\n'.join(config_list))

    fn_suffix = f"N({N})_alpha({alpha:0.4f})_seed({seed})"

    fn_test = data_path + '/lindbladian_evals_' + fn_suffix + '.txt'

    if not os.path.isfile(fn_test):
        print('file does not exist')
        if segment == 'short':
            os.system('sbatch run_short.sh ' + data_path)
        elif segment == 'medium':
            os.system('sbatch run_medium.sh ' + data_path)
        elif segment == 'long':
            os.system('sbatch run_long.sh ' + data_path)
