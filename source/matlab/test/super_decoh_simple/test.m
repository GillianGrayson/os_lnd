clear all;

reshuffle_type = 1; % 0 - original, 1 - mine
G_type = 0;
N = 10; % system size
p = 0.1; % parameter for partial decoherence 0<=p<=1
seed = 1;

check_f_basis = 1;

cpp_path = 'E:/Work/os_lnd/source/cpp/os_lnd/os_lnd';
suffix = sprintf('reshuffle(%d)_G(%d)_N(%d)_p(%0.4f)_seed(%d)', reshuffle_type, G_type, N, p, seed);

N2 = N^2;
M = N^2; % auxiliary size
imag1 = sqrt(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = zeros(M);
fn_cpp = sprintf('%s/G_mtx_%s.txt', cpp_path, suffix);
cpp_data = importdata(fn_cpp);
for st_id_1 = 1:M
    for st_id_2 = 1:M
        str_id = (st_id_1 - 1) * M + st_id_2;
        str = string(cpp_data(str_id));
        data = sscanf(str, '(%e,%e)', 2);
        G(st_id_1, st_id_2) = data(1) + 1i * data(2);
    end
end
G = transpose(G); % Eigen is column-major by default.Here we have row-major dense matrix

rho = Decoh(G,M,p); % decoherence of the random state

if reshuffle_type == 1
    F = reshuffle_1(rho,N); % reshuffling procedure
else
    F = reshuffle_0(rho,N); % reshuffling procedure
end
    
y1 = zeros(M,1);
y2 = eye(N);
for s1 = 1:N
    for s2 = 1:N
        kk=s2+(N-1)*s1;
        y1(kk)=y2(s1,s2);
    end
end

FT=ctranspose(F);

S=FT*y1;

S2 = zeros(N,N);
for s1 = 1:N
    for s2 = 1:N
        kk=s2+(N-1)*s1;
        S2(s1,s2)=S(kk);
    end
end

trace(S2)

P = F - (1/2)*(kron(S2,eye(N)) + kron(eye(N),transpose(S2)));

L = zeros(N2);
fn_cpp = sprintf('%s/lindbladian_mtx_%s.txt', cpp_path, suffix);
cpp_data = importdata(fn_cpp);
for st_id_1 = 1:N2
    for st_id_2 = 1:N2
        str_id = (st_id_1 - 1) * N2 + st_id_2;
        str = string(cpp_data(str_id));
        data = sscanf(str, '(%e,%e)', 2);
        L(st_id_1, st_id_2) = data(1) + 1i * data(2);
    end
end
L = transpose(L); % Eigen is column-major by default. Here we have row-major dense matrix

lindbladian_diff = norm(L - P)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Lindbladian evals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Evec,D] = eig(P); % eigenvalues of Lindblad superoperator
Eval = diag(D);
[C, min_eval_id] = min(abs(Eval));
Eval_sorted = sort(Eval, 'ComparisonMethod', 'abs');

evals = zeros(N2, 1);
fn_cpp = sprintf('%s/lindbladian_evals_%s.txt', cpp_path, suffix);
cpp_data = importdata(fn_cpp);
for str_id = 1:N2
        str = string(cpp_data(str_id));
        data = sscanf(str, '(%e,%e)', 2);
        evals(str_id) = data(1) + 1i * data(2);
end
evals_sorted = sort(evals, 'ComparisonMethod', 'abs');

evals_diff = norm(abs(evals_sorted) - abs(Eval_sorted))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Rho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho_ml = zeros(N,N);

for i=1:N
    for j=1:N
        rho_ml(i,j) = Evec((i-1)*N+j,min_eval_id);
    end
end

rho_ml = rho_ml/trace(rho_ml);

rho = zeros(N);
fn_cpp = sprintf('%s/rho_mtx_%s.txt', cpp_path, suffix);
cpp_data = importdata(fn_cpp);
for st_id_1 = 1:N
    for st_id_2 = 1:N
        str_id = (st_id_1 - 1) * N + st_id_2;
        str = string(cpp_data(str_id));
        data = sscanf(str, '(%e,%e)', 2);
        rho(st_id_1, st_id_2) = data(1) + 1i * data(2);
    end
end
rho = transpose(rho); % Eigen is column-major by default.Here we have row-major dense matrix

rho_diff = norm(rho_ml - rho)