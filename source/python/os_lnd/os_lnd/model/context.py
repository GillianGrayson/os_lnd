from os_lnd.model.strategy.setup.mbl import MBLSetupStrategy

class Context:
    def __init__(self, model):
        system = model.ini['global']['system']
        if system == 'mbl':
            self.setup_strategy = MBLSetupStrategy()
        else:
            raise ValueError(f'Unsupported system: {system}')

    def run(self, model):
        self.setup_strategy.setup_model(model)
        model.calc_zero_eigen_vector()