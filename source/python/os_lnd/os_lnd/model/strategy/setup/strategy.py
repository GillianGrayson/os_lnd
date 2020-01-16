import abc
from scipy import sparse
import numpy as np

class SetupStrategy(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def setup_model(self, model):
        pass


    def calc_lindbladian(self, model):

        eye = sparse.eye(model.sys_size, model.sys_size, dtype=np.complex, format='csr')

        model.lindbladian = -1.0j * (
                    sparse.kron(eye, model.hamiltonian) -
                    sparse.kron(model.hamiltonian.transpose(copy=True), eye)
        )

        for diss_id, diss in enumerate(model.dissipators):

            tmp_1 = diss.getH().transpose(copy=True)
            tmp_2 = diss.getH() * diss
            tmp_3 = tmp_2.transpose(copy=True)

            model.lindbladian += 0.5 * model.diss_gammas[diss_id] * (
                        2.0 * sparse.kron(eye, diss) * sparse.kron(tmp_1, eye) -
                        sparse.kron(tmp_3, eye) -
                        sparse.kron(eye, tmp_2)
            )


