clear all;

Nc                  = 8;
seed                = 10;
dissipator_type     = 1;
alpha               = 0;
g                   = 0.1;
W                   = 3;
U                   = 1;
J                   = 1;

tic

cpp_path = 'E:/Work/os_lnd/source/cpp/os_lnd/os_lnd';

Np = Nc/2;
Ns = nchoosek(Nc,Np);
num_states = precalc_states(Nc, Np, 0);

file_name_suffix = sprintf('ns(%d)_seed(%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)', ...
    Nc, ...
    seed, ...
    dissipator_type, ...
    alpha, ...
    g, ...
    W, ...
    U, ...
    J);

H  = zeros(Ns); % Halimtonian
Hd = zeros(Ns); % Disorder
Hi = zeros(Ns); % Interaction
Hh = zeros(Ns); % Hopping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Disorder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_name_cpp = sprintf('%s/energies_%s.txt', cpp_path, file_name_suffix);
E = importdata(file_name_cpp);

for state_id_1 = 1:Ns
    Hd(state_id_1, state_id_1) = ( dec2bin(idtox(state_id_1), Nc) == '1' ) * E;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Interaction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for state_id_1 = 1:Ns
    shifted = bitshift(idtox(state_id_1), 1);
    Hi(state_id_1, state_id_1) = sum( dec2bin(bitand(idtox(state_id_1), shifted)) == '1' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Hopping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for state_id_1 = 1:Ns
    for state_id_2 = 1:Ns
        Hh(state_id_1, state_id_2) = is_adjacent(idtox(state_id_1), idtox(state_id_2));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = -J*Hh + U*Hi + 1.0*W*Hd;

cpp_hamitlonian = zeros(Ns, Ns);
file_name_cpp = sprintf('%s/hamiltonian_mtx_%s.txt', cpp_path, file_name_suffix);
cpp_hamitlonian_data = importdata(file_name_cpp);
for str_id = 1 : size(cpp_hamitlonian_data, 1)
    str = string(cpp_hamitlonian_data(str_id));
    data = sscanf(str, '%d\t%d\t(%e,%e)',4);
    cpp_hamitlonian(data(1) + 1, data(2) + 1) = data(3) + 1i * data(4);
end
hamiltonian_check = max(max(abs(H - cpp_hamitlonian)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Creating supermatrix of right part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

super_rp_matrix = -sqrt(-1) * ( kron(eye(Ns), H) - kron(transpose(H), eye(Ns)) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Dissipator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(dissipator_type == 1)
    
    for dissipator_id = 1:Nc-1
        
        dissipator = zeros(Ns);
        
        for state_id_1 = 1:Ns
            dissipator(state_id_1, state_id_1) = bitget(idtox(state_id_1), dissipator_id) - bitget(idtox(state_id_1), dissipator_id+1);
            for state_id_2 = 1:Ns
                if(is_adjacent(idtox(state_id_1), idtox(state_id_2)))
                    hopping_ids = 1 + Nc - find( dec2bin(bitxor(idtox(state_id_1), idtox(state_id_2)),Nc) == '1' );
                    if(hopping_ids(2) == dissipator_id)
                        if(bitget(idtox(state_id_1), dissipator_id))
                            dissipator(state_id_2,state_id_1) = exp(sqrt(-1)*alpha);
                            dissipator(state_id_1,state_id_2) = -exp(-sqrt(-1)*alpha);
                        else
                            dissipator(state_id_2,state_id_1) = -exp(-sqrt(-1)*alpha);
                            dissipator(state_id_1,state_id_2) = exp(sqrt(-1)*alpha);
                        end
                    end
                end
            end
        end
        
        super_rp_matrix = super_rp_matrix + ...
            g * 0.5 * (2.0 * kron(eye(Ns),dissipator) * kron(transpose(dissipator'),eye(Ns)) - ...
            kron(transpose(dissipator'*dissipator), eye(Ns)) - ...
            kron(eye(Ns), dissipator'*dissipator));
        
    end
    
elseif(dissipator_type == 0)
    
    for dissipator_id = 1:Nc
        
        dissipator=zeros(Ns);
        for state_id_1=1:Ns
            dissipator(state_id_1, state_id_1) = bitget(idtox(state_id_1), dissipator_id);
        end

        super_rp_matrix = super_rp_matrix + ...
            g * 0.5 * (2.0 * kron(eye(Ns),dissipator) * kron(transpose(dissipator'),eye(Ns)) - ...
            kron(transpose(dissipator'*dissipator), eye(Ns)) - ...
            kron(eye(Ns), dissipator'*dissipator));
    end
    
end

cpp_lind = zeros(Ns, Ns);
file_name_cpp = sprintf('%s/lindbladian_mtx_%s.txt', cpp_path, file_name_suffix);
cpp_lind_data = importdata(file_name_cpp);
for str_id = 1 : size(cpp_lind_data, 1)
    str = string(cpp_lind_data(str_id));
    data = sscanf(str, '%d\t%d\t(%e,%e)',4);
    cpp_lind(data(1) + 1, data(2) + 1) = data(3) + 1i * data(4);
end
lind_check = max(max(abs(super_rp_matrix - cpp_lind)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Zero eigen vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sparse_super_rp_matrix = sparse(super_rp_matrix);
super_rp_matrix = 0;

[zero_evec, zero_eval] = eigs(sparse_super_rp_matrix, 1, 'sm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Rho in direct basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho = zeros(Ns, Ns);
for state_id_1 = 1:Ns
    rho(:, state_id_1) = zero_evec(1+(state_id_1-1)*(Ns) : state_id_1*(Ns));
end
rho = rho / trace(rho);

cpp_rho = zeros(Ns, Ns);
file_name_cpp = sprintf('%s/rho_mtx_%s.txt', cpp_path, file_name_suffix);
cpp_rho_data = importdata(file_name_cpp);
for st_id_1 = 1:Ns
    for st_id_2 = 1:Ns
        str_id = (st_id_1 - 1) * Ns + st_id_2;
        str = string(cpp_rho_data(str_id));
        data = sscanf(str, '(%e,%e)',4);
        cpp_rho(st_id_1, st_id_2) = data(1) + 1i * data(2);
    end
end
rho_check = max(max(abs(rho - cpp_rho)))
rho_check_trans = max(max(abs(rho - transpose(cpp_rho))))
