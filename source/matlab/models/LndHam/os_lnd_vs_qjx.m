clear all;

N = 10;
alpha = 0.5;
seed = 1;

qjx_T = 1;

os_lnd_path = 'E:/Work/os_lnd/source/cpp/os_lnd/os_lnd';
qjx_path = 'E:/Work/os_d/source/cpp/QJX/QJX';

os_lnd_suffix = sprintf('N(%d)_alpha(%0.4f)_seed(%d)', ...
    N, ...
    alpha, ...
    seed);

qjx_suffix = sprintf('setup(4_1_0)_rnd(1_1000000)_N(%d)_rnd(%d)_alpha(%0.4f)_T(%0.4f)', ...
    N, ...
    seed, ...
    alpha, ...
    qjx_T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = sprintf('%s/hamiltonian_mtx_%s.txt', os_lnd_path, os_lnd_suffix);
os_lnd_data = importdata(fn);
os_lnd_H = zeros(N);
for st_id_1 = 1:N
    for st_id_2 = 1:N
        str_id = (st_id_1 - 1) * N + st_id_2;
        str = string(os_lnd_data(str_id));
        data = sscanf(str, '(%e,%e)', 2);
        os_lnd_H(st_id_1, st_id_2) = data(1) + 1i * data(2);
    end
end
os_lnd_H = transpose(os_lnd_H); % Eigen is column-major by default.Here we have row-major dense matrix
os_lnd_H = alpha / sqrt(N) * os_lnd_H;

fn = sprintf('%s/hamiltonian_%s.txt', qjx_path, qjx_suffix);
qjx_data = importdata(fn);
qjx_H = zeros(N);
for s_id_1 = 1:N
    for s_id_2 = 1:N
        index = (s_id_1 - 1) * N + s_id_2;
        qjx_H(s_id_1, s_id_2) = qjx_data(index, 1) + 1i * qjx_data(index, 2);
    end
end

H_check = norm(os_lnd_H - qjx_H)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 Rho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = sprintf('%s/rho_mtx_%s.txt', os_lnd_path, os_lnd_suffix);
os_lnd_data = importdata(fn);
os_lnd_rho = zeros(N, N);
for st_id_1 = 1:N
    for st_id_2 = 1:N
        str_id = (st_id_1 - 1) * N + st_id_2;
        str = string(os_lnd_data(str_id));
        data = sscanf(str, '(%e,%e)',4);
        os_lnd_rho(st_id_1, st_id_2) = data(1) + 1i * data(2);
    end
end
os_lnd_diag_rho = abs(diag(os_lnd_rho));

fn = sprintf('%s/adr_avg_%s.txt', qjx_path, qjx_suffix);
qjx_data = importdata(fn);

rho_diag_diff = norm(os_lnd_diag_rho - qjx_data)
