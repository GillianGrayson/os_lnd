from scipy.special import comb
from os_lnd.model.strategy.setup.strategy import SetupStrategy
from os_lnd.infrastructure.load import load_sp_mtx
from os_lnd.infrastructure.file_system import get_input_path


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
        model.diss_gammas = []
        for diss_id in range(0, model.num_dissipators):
            fn = get_input_path(model.params_path) + f'/diss_{diss_id}_mtx' + model.suffix
            model.dissipators.append(load_sp_mtx(fn, model.sys_size))
            model.diss_gammas.append(diss_gamma)

        self.calc_lindbladian(model)



