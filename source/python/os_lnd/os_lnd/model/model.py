class Model:
    def __init__(self, ini):
        self.ini = ini

        self.setup_strategy = None
        self.params_path = None
        self.suffix = None
        self.sys_size = None
        self.num_dissipators = None
        self.hamiltonian = None
        self.dissipators = None
        self.lindbladian = None
