from scipy.sparse.linalg import eigs

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

    def calc_zero_eigen_vector(self):
        vals, vecs = eigs(self.lindbladian, k=1, which='SM')

        ololo = 1