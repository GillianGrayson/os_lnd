clear all;

N = 5;

tau = 1;
k = -1;

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
%                  All middle dissipators (N-1, N)
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

% =====================================================================
%      Check the first dissipator which are constructed this way
% =====================================================================
first_diss = zeros(num_states);
for state_id = 1:num_states
    state_int = state_id - 1;
    state_bits = de2bi(state_int, N);
    state_bits = flip(state_bits);
    
    if (state_bits(1) == state_bits(2))
        if state_bits(1) == 0
            first_diss(state_id, state_id) = tau;
        else
            first_diss(state_id, state_id) = k;
        end
    else
        inv_bits = state_bits;
        inv_bits(1) = state_bits(2);
        inv_bits(2) = state_bits(1);
        inv_bits = flip(inv_bits);
        row_id = bi2de(inv_bits) + 1;
        first_diss(row_id, state_id) = 1;
    end
end
fprintf('first dissipator difference = %0.4e\n', norm(dissipators{1} - first_diss));

% =====================================================================
%                            Lindbladian
% =====================================================================
Lind = zeros(num_states^2);
for diss_id = 1:size(dissipators, 1)
    tmp = 0.5 * (2 * kron(eye(num_states), dissipators{diss_id}) * ...
        kron(transpose(dissipators{diss_id}'), eye(num_states)) - ...
        kron(transpose(dissipators{diss_id}' * dissipators{diss_id}), eye(num_states)) - ...
        kron(eye(num_states), dissipators{diss_id}' * dissipators{diss_id}));
    Lind = Lind + tmp;
end

% =====================================================================
%                            Evals routines
% =====================================================================
evals = eig(Lind);
