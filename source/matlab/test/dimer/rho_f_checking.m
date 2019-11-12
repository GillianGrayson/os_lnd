clear all;

global N J E0 A0 w U g HE HU HJ A phi0 Ps Pd

imag1=sqrt(-1);

N=10; % number of particles, system size = N+1
Ns = N + 1;

E0 = 0.; % bias
U = 0.5; % on-site interaction
J = -1; % hopping constant

g=0.1/N; % gamma;

A0 = -3.4; %driving amplitude
w = 1; % driving frequency
phi0 = 0; % initial phase

cpp_path = 'E:/Work/os_lnd/source/cpp/os_lnd/os_lnd';

file_name_suffix = sprintf('np(%d)_diss(1_0.1000)_prm(%0.4f_%0.4f_%0.4f)_drv(1_%0.4f_%0.4f_%0.4f)', ...
    N, ...
    E0, ...
    U, ...
    -J, ...
    -A0, ...
    w, ...
    phi0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HE=zeros(N+1);
HJ=zeros(N+1);
HU=zeros(N+1);

HE=diag(2*[0:N]-N); % diagonal
HU=diag((2*[0:N]-N).^2); % diagonal
for i=1:N
    HJ(i+1,i)=sqrt((N-i+1)*(i));
    HJ(i,i+1)=sqrt((i)*(N-i+1));
end

HE=HE-trace(HE)/(N+1)*eye(N+1);
HU=HU-trace(HU)/(N+1)*eye(N+1);
HJ=HJ-trace(HJ)/(N+1)*eye(N+1);

H=(E0*HE+U/N*HU+J*HJ);

cpp_hamitlonian = zeros(Ns, Ns);
file_name_cpp = sprintf('%s/hamiltonian_mtx_%s.txt', cpp_path, file_name_suffix);
cpp_hamitlonian_data = importdata(file_name_cpp);
for str_id = 1 : size(cpp_hamitlonian_data, 1)
    str = string(cpp_hamitlonian_data(str_id));
    data = sscanf(str, '%d\t%d\t(%e,%e)',4);
    cpp_hamitlonian(data(1) + 1, data(2) + 1) = data(3) + 1i * data(4);
end
hamiltonian_check = max(max(abs(H - cpp_hamitlonian)))

cpp_hamitlonian = zeros(Ns, Ns);
file_name_cpp = sprintf('%s/hamiltonian_drv_mtx_%s.txt', cpp_path, file_name_suffix);
cpp_hamitlonian_data = importdata(file_name_cpp);
for str_id = 1 : size(cpp_hamitlonian_data, 1)
    str = string(cpp_hamitlonian_data(str_id));
    data = sscanf(str, '%d\t%d\t(%e,%e)',4);
    cpp_hamitlonian(data(1) + 1, data(2) + 1) = data(3) + 1i * data(4);
end
hamiltonian_drv_check = max(max(abs(HE - cpp_hamitlonian)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dissipator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A1=diag(2*[0:N]-N);
A2=zeros(N+1);
for i=1:N
    A2(i+1,i)=-sqrt((N-i+1)*(i));
    A2(i,i+1)=sqrt((i)*(N-i+1));
end
A2=imag1*A2;

A=A1-sqrt(-1)*A2;

cpp_mtx = zeros(Ns, Ns);
file_name_cpp = sprintf('%s/diss_type_0_mtx_%s.txt', cpp_path, file_name_suffix);
cpp_mtx_data = importdata(file_name_cpp);
for str_id = 1 : size(cpp_mtx_data, 1)
    str = string(cpp_mtx_data(str_id));
    data = sscanf(str, '%d\t%d\t(%e,%e)',4);
    cpp_mtx(data(1) + 1, data(2) + 1) = data(3) + 1i * data(4);
end
dissipator_check = max(max(abs(A - cpp_mtx)))

P=zeros((N+1)^2);
P=-sqrt(-1)*(kron(eye(N+1),H)-kron(transpose(H),eye(N+1)))+...
    g/2*(2*kron(eye(N+1),A)*kron(transpose(A'),eye(N+1))-...
    kron(transpose(A'*A),eye(N+1))-kron(eye(N+1),A'*A));
Ps=sparse(P);
P=0;

cpp_mtx = zeros(Ns^2);
file_name_cpp = sprintf('%s/lindbladian_mtx_%s.txt', cpp_path, file_name_suffix);
cpp_mtx_data = importdata(file_name_cpp);
for str_id = 1 : size(cpp_mtx_data, 1)
    str = string(cpp_mtx_data(str_id));
    data = sscanf(str, '%d\t%d\t(%e,%e)',4);
    cpp_mtx(data(1) + 1, data(2) + 1) = data(3) + 1i * data(4);
end
lindbladian_mtx_check = max(max(abs(Ps - cpp_mtx)))

Pd=sparse(-sqrt(-1)*(kron(eye(N+1),HE)-kron(transpose(HE),eye(N+1))));

cpp_mtx = zeros(Ns^2);
file_name_cpp = sprintf('%s/lindbladian_drv_mtx_%s.txt', cpp_path, file_name_suffix);
cpp_mtx_data = importdata(file_name_cpp);
for str_id = 1 : size(cpp_mtx_data, 1)
    str = string(cpp_mtx_data(str_id));
    data = sscanf(str, '%d\t%d\t(%e,%e)',4);
    cpp_mtx(data(1) + 1, data(2) + 1) = data(3) + 1i * data(4);
end
lindbladian_drv_mtx_check = max(max(abs(Pd - cpp_mtx)))


rho_mtx_before = zeros(N+1);
rho_mtx_before(1,1) = 1;
rho_before = zeros((N+1)^2, 1);
for st_1 = 1:N+1
    for st_2 = 1:N+1
        index = (st_1-1) * (N+1) + st_2;
        rho_before(index) = rho_mtx_before(st_1, st_2);
    end
end

t0 = 0.0;
t1 = 10.0 * 2*pi/w;

options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[t,Y]=ode45('sysLind',[t0:pi/w:t1], rho_before, options);
rho_after = transpose(Y(length(t),:));

rho_after_mtx = zeros(N+1);
for st_1 = 1:N+1
    for st_2 = 1:N+1
        index = (st_1-1) * (N+1) + st_2;
        rho_after_mtx(st_1, st_2) = rho_after(index);
    end
end

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
rho_check = max(max(abs(rho_after_mtx - cpp_rho)))
rho_check_trans = max(max(abs(rho_after_mtx - transpose(cpp_rho))))


