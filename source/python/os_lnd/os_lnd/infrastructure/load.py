import numpy as np
from tqdm import tqdm
from scipy.sparse import csr_matrix
import re


def num_lines_in_file(fn):
    with open(fn) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def load_sp_mtx(fn, size):

    num_lines = num_lines_in_file(fn)

    rows = np.zeros(num_lines, dtype=np.int)
    cols = np.zeros(num_lines, dtype=np.int)
    data = np.zeros(num_lines, dtype=np.complex)

    f = open(fn)
    for line_id, line in tqdm(enumerate(f), mininterval=60.0, desc='load_sp_mtx'):
        m = re.match(r'(?P<row>.*)\t(?P<col>.*)\t\((?P<real>.*),(?P<imag>.*)\)', line)
        rows[line_id] = int(m.group('row'))
        cols[line_id] = int(m.group('col'))
        data[line_id] = float(m.group('real')) + float(m.group('imag')) * 1j
    f.close()

    mtx = csr_matrix((data, (rows, cols)), shape=(size, size))

    return mtx


def load_mtx(fn, size):
    num_lines = num_lines_in_file(fn)

    data = np.zeros(num_lines, dtype=np.complex)
    f = open(fn)
    for line_id, line in tqdm(enumerate(f), mininterval=60.0, desc='load_mtx'):
        m = re.match(r'\((?P<real>.*),(?P<imag>.*)\)', line)
        data[line_id] = float(m.group('real')) + float(m.group('imag')) * 1j
    f.close()

    mtx = data.reshape((size,size))

    return mtx