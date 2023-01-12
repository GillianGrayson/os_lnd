clear all;

Nc = 6;
seed = 1;
diss_type = 1;
diss_phase = 0;
diss_gamma = 0.1;
W = 1;
U = 1;
J = 1;

os_lnd_path = 'E:/Work/os_lnd/source/cpp/os_lnd/os_lnd';
qjx_path = 'E:/Work/os_d/source/cpp/QJX/QJX';

Np = Nc/2;
Ns = nchoosek(Nc,Np);

os_lnd_suffix = sprintf('ns(%d)_seed(%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)', ...
    Nc, ...
    seed, ...
    diss_type, ...
    diss_phase, ...
    diss_gamma, ...
    W, ...
    U, ...
    J);

qjx_suffix = sprintf('setup(3_1_0)_rnd(1_1000000)_Nc(%d)_rnd(%d_1000000)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(0_0)', ...
    Nc, ...
    seed, ...
    diss_type, ...
    diss_phase, ...
    diss_gamma, ...
    W, ...
    U, ...
    J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = sprintf('%s/hamiltonian_mtx_%s.txt', os_lnd_path, os_lnd_suffix);
os_lnd_data = importdata(fn);
os_lnd_H = zeros(Ns, Ns);
for str_id = 1 : size(os_lnd_data, 1)
    str = string(os_lnd_data(str_id));
    data = sscanf(str, '%d\t%d\t(%e,%e)',4);
    os_lnd_H(data(1) + 1, data(2) + 1) = data(3) + 1i * data(4);
end

fn = sprintf('%s/hamiltonian_%s.txt', qjx_path, qjx_suffix);
qjx_data = importdata(fn);
qjx_H = zeros(Ns, Ns);
for s_id_1 = 1:Ns
    for s_id_2 = 1:Ns
        index = (s_id_1 - 1) * Ns + s_id_2;
        qjx_H(s_id_1, s_id_2) = qjx_data(index);
    end
end

H_check = norm(os_lnd_H - qjx_H)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Dissipators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for diss_id = 1:Nc-1
    
    fn = sprintf('%s/diss_%d_mtx_%s.txt', os_lnd_path, diss_id - 1, os_lnd_suffix);
    os_lnd_data = importdata(fn);
    os_lnd_diss = zeros(Ns, Ns);
    for str_id = 1 : size(os_lnd_data, 1)
        str = string(os_lnd_data(str_id));
        data = sscanf(str, '%d\t%d\t(%e,%e)',4);
        os_lnd_diss(data(1) + 1, data(2) + 1) = data(3) + 1i * data(4);
    end
    
    fn = sprintf('%s/dissipator_%d_%s.txt', qjx_path, diss_id - 1, qjx_suffix);
    qjx_data = importdata(fn);
    qjx_diss = zeros(Ns, Ns);
    for s_id_1 = 1:Ns
        for s_id_2 = 1:Ns
            index = (s_id_1 - 1) * Ns + s_id_2;
            qjx_diss(s_id_1, s_id_2) = qjx_data(index);
        end
    end
    
    diss_check = norm(os_lnd_diss - qjx_diss)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 Rho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = sprintf('%s/rho_mtx_%s.txt', os_lnd_path, os_lnd_suffix);
os_lnd_data = importdata(fn);
os_lnd_rho = zeros(Ns, Ns);
for st_id_1 = 1:Ns
    for st_id_2 = 1:Ns
        str_id = (st_id_1 - 1) * Ns + st_id_2;
        str = string(os_lnd_data(str_id));
        data = sscanf(str, '(%e,%e)',4);
        os_lnd_rho(st_id_1, st_id_2) = data(1) + 1i * data(2);
    end
end
os_lnd_diag_rho = abs(diag(os_lnd_rho));

fn = sprintf('%s/adr_avg_%s.txt', qjx_path, qjx_suffix);
qjx_data = importdata(fn);

rho_diag_diff = norm(os_lnd_diag_rho - qjx_data)
