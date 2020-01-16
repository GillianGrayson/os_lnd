import abc


class SetupStrategy(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def setup_model(self, model):
        pass



