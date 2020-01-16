from scipy.sparse.linalg import eigs
import numpy as np
from os_lnd.infrastructure.load import load_mtx
from os_lnd.infrastructure.file_system import get_input_path
from scipy import linalg as sla
from numpy import linalg


class Model:
    def __init__(self, ini):
        self.ini = ini

        self.setup_strategy = None
        self.params_path = None
        self.suffix = None
        self.sys_size = None
        self.num_dissipators = None
        self.diss_gammas = None
        self.hamiltonian = None
        self.dissipators = None
        self.lindbladian = None

        self.rho = None
        self.lindbladian_evals = None

    def calc_zero_eigen_vector(self):
        vals, vecs = eigs(self.lindbladian, k=1, which='SM')

        rho = vecs.reshape((self.sys_size, self.sys_size))
        trace = np.trace(rho)
        rho = rho / trace
        trace = np.trace(rho)
        print(f'trace = {trace}')

        self.rho = rho

    def calc_lindbladian_evals(self):
        aaa = eigs(self.lindbladian, k=self.sys_size * self.sys_size - 2, return_eigenvectors=False, which='LR')
        bbb = eigs(self.lindbladian, k=2, return_eigenvectors=False, which='SR')
        print(f'rho_difs')

    def check_rho(self):
        fn = get_input_path(self.params_path) + '/rho_mtx' + self.suffix
        rho_original = load_mtx(fn, self.sys_size)

        rho_diff = self.rho - rho_original

        print(f'rho_diff = {linalg.norm(rho_diff)}')