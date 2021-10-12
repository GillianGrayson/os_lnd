clear all;

N = 5;

tau = 1;
k = -1;

qjx_seed = 1;
qjx_T = 1;
qjx_path = 'E:/Work/os_d/source/cpp/QJX/QJX';
qjx_suffix = sprintf('setup(5_8_0)_rnd(1_1000000)_N(%d)_seed(%d)_tau(%d)_k(%d)_T(%0.4f)_lpn(-1_-3.0000_-3.0000_-3.0000)', ...
    N, ...
    qjx_seed, ...
    tau, ...
    k, ...
    qjx_T);


con = zeros(4);
con(1, 1) = tau;
con(4, 4) = k;
con(2, 3) = 1;
con(3, 2) = 1;

num_states = 2^N;

dissipators = cell(N, 1);

% =====================================================================
%                       Left dissipator (1,2)
% =====================================================================
left_diss = con;
for spin_id = 3:N
    left_diss = kron(left_diss, eye(2));
end
dissipators{1} = left_diss;

% =====================================================================
%                      Right dissipator (N-1, N)
% =====================================================================
right_diss = eye(2);
for spin_id = 2:N-2
    right_diss = kron(right_diss, eye(2));
end
right_diss = kron(right_diss, con);
dissipators{N-1} = right_diss;

% =====================================================================
%                  All middle dissipators
% =====================================================================
for diss_id = 2:N-2
    
    curr_diss = eye(2);
    
    left_s_id = 2;
    left_f_id = diss_id - 1;
    for spin_id = left_s_id : left_f_id
        curr_diss = kron(curr_diss, eye(2));
    end
    
    curr_diss = kron(curr_diss, con);
    
    right_s_id = diss_id + 2;
    right_f_id = N;
    for spin_id = right_s_id : right_f_id
        curr_diss = kron(curr_diss, eye(2));
    end
    
    dissipators{diss_id} = curr_diss;
end

% =====================================================================
%                       Border dissipator (N, 1)
% =====================================================================
border_diss = zeros(num_states);
for state_id = 1:num_states
    state_int = state_id - 1;
    state_bits = de2bi(state_int, N);
    state_bits = flip(state_bits);
    
    if (state_bits(1) == state_bits(end))
        if state_bits(1) == 0
            border_diss(state_id, state_id) = tau;
        else
            border_diss(state_id, state_id) = k;
        end
    else
        inv_bits = state_bits;
        inv_bits(1) = state_bits(end);
        inv_bits(end) = state_bits(1);
        inv_bits = flip(inv_bits);
        row_id = bi2de(inv_bits) + 1;
        border_diss(row_id, state_id) = 1;
    end
end
dissipators{N} = border_diss;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = sprintf('%s/hamiltonian_%s.txt', qjx_path, qjx_suffix);
qjx_data = importdata(fn);
qjx_H = zeros(num_states);
for s_id_1 = 1:num_states
    for s_id_2 = 1:num_states
        index = (s_id_1 - 1) * num_states + s_id_2;
        qjx_H(s_id_1, s_id_2) = qjx_data(index, 1) + 1i * qjx_data(index, 2);
    end
end
H_check = norm(qjx_H)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Dissipators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for diss_id = 1:size(dissipators, 1)
    fn = sprintf('%s/dissipator_%d_%s.txt', qjx_path, diss_id - 1, qjx_suffix);
    qjx_data = importdata(fn);
    qjx_d = zeros(num_states);
    for s_id_1 = 1:num_states
        for s_id_2 = 1:num_states
            index = (s_id_1 - 1) * num_states + s_id_2;
            qjx_d(s_id_1, s_id_2) = qjx_data(index, 1) + 1i * qjx_data(index, 2);
        end
    end
    d_check = norm(qjx_d - dissipators{diss_id})
end

