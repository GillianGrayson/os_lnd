from scipy.special import comb
from os_lnd.model.strategy.setup.strategy import SetupStrategy
from os_lnd.infrastructure.load import load_sp_mtx
from os_lnd.infrastructure.file_system import get_input_path
from scipy import sparse
from scipy.sparse.linalg import norm
import numpy as np
import copy


class MBLSetupStrategy(SetupStrategy):
    def setup_model(self, model):

        name_precision = model.ini['global']['name_precision']

        num_spins = int(model.ini['mbl']['num_spins'])
        seed = int(model.ini['mbl']['seed'])

        diss_type = int(model.ini['mbl']['diss_type'])
        diss_phase = float(model.ini['mbl']['diss_phase'])
        diss_gamma = float(model.ini['mbl']['diss_gamma'])

        W = float(model.ini['mbl']['W'])
        U = float(model.ini['mbl']['U'])
        J = float(model.ini['mbl']['J'])

        model.params_path = ''
        precision = f'.{name_precision}f'
        model.suffix = f'_ns({num_spins})' + \
                       f'_seed({seed})' + \
                       f'_diss({diss_type}_{diss_phase:{precision}}_{diss_gamma:{precision}})' + \
                       f'_prm({W:{precision}}_{U:{precision}}_{J:{precision}})' + \
                       '.txt'

        model.sys_size = comb(num_spins, num_spins / 2, True)

        if diss_type == 0:
            model.num_dissipators = num_spins
        elif diss_type == 1:
            model.num_dissipators = num_spins - 1
        else:
            raise ValueError(f'Unsupported dissipator type')

        fn = get_input_path(model.params_path) + '/hamiltonian_mtx' + model.suffix
        model.hamiltonian = load_sp_mtx(fn, model.sys_size)

        model.dissipators = []
        for diss_id in range(0, model.num_dissipators):
            fn = get_input_path(model.params_path) + f'/diss_{diss_id}_mtx' + model.suffix
            model.dissipators.append(load_sp_mtx(fn, model.sys_size))

        eye = sparse.eye(model.sys_size, model.sys_size, dtype=np.complex, format='csr')

        model.lindbladian = -1.0j * (sparse.kron(eye, model.hamiltonian) - sparse.kron(model.hamiltonian.transpose(copy=True), eye))

        for diss in model.dissipators:

            tmp_1 = diss.getH().transpose(copy=True)
            tmp_2 = diss.getH() * diss
            tmp_3 = tmp_2.transpose(copy=True)

            model.lindbladian += 0.5 * diss_gamma * (2.0 * sparse.kron(eye, diss) * sparse.kron(tmp_1, eye) - sparse.kron(tmp_3, eye) - sparse.kron(eye, tmp_2))

        fn = get_input_path(model.params_path) + '/lindbladian_mtx' + model.suffix
        lindbladian_origin = load_sp_mtx(fn, model.sys_size * model.sys_size)

        lindbladian_diff = lindbladian_origin - model.lindbladian

        n = norm(lindbladian_diff)

        ololo = 1



