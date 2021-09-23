clear all;


N = 10; % system size
alpha = 0.5;
seed = 1;

check_f_basis = 1;

cpp_path = 'E:/Work/os_lnd/source/cpp/os_lnd/os_lnd';
suffix = sprintf('N(%d)_alpha(%0.4f)_seed(%d)', N, alpha, seed);

M = N^2-1; % auxiliary size
N2 = N^2;

%%% Step 1
%%% NB: here we complete the F-basis by a normalized identity matrix as
%%% F{1} - just for the future, but will not actually use it.
k=1;
F{1}=sparse(eye(N))/sqrt(N);
for i=1:N
    for j=i+1:N
        k=k+1;
        F{k}=sparse([i j],[j i],[1 1]/sqrt(2),N,N);
        k=k+1;
        F{k}=sparse([i j],[j i],-1i*[1 -1]/sqrt(2),N,N);
    end
end

for i=1:N-1
    k=k+1;
    temp=zeros(1,i+1);
    temp(1:i)=ones(1,i);
    temp(i+1)=-i;
    F{k}=sparse([1:i+1],[1:i+1],temp/sqrt(i*(i+1)),N,N);
end


if (check_f_basis)
    diffs = zeros(size(F, 2), 1);
    for fb_id = 1:size(F, 2)
        test_mtx = zeros(N, N);
        fn_cpp = sprintf('%s/f_basis_%d_mtx_%s.txt', cpp_path, fb_id - 1, suffix);
        cpp_data = importdata(fn_cpp);
        for str_id = 1 : size(cpp_data, 1)
            str = string(cpp_data(str_id));
            data = sscanf(str, '%d\t%d\t(%e,%e)', 4);
            test_mtx(data(1) + 1, data(2) + 1) = data(3) + 1i * data(4);
        end
        curr_f = full(F{fb_id});
        diffs(fb_id) = norm(curr_f - test_mtx);
    end
    max_fb_diff = max(diffs)
end


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ham = zeros(N);
fn_cpp = sprintf('%s/hamiltonian_mtx_%s.txt', cpp_path, suffix);
cpp_data = importdata(fn_cpp);
for st_id_1 = 1:N
    for st_id_2 = 1:N
        str_id = (st_id_1 - 1) * N + st_id_2;
        str = string(cpp_data(str_id));
        data = sscanf(str, '(%e,%e)', 2);
        Ham(st_id_1, st_id_2) = data(1) + 1i * data(2);
    end
end
Ham = transpose(Ham); % Eigen is column-major by default.Here we have row-major dense matrix


[Evec,D] = eig(G);
for k1=1:M
    DISS{k1}=0;
    for k2=1:M
        DISS{k1}=DISS{k1} + Evec(k2,k1) * F{k2+1};
    end
end

%%% Step 3: calculate H-matrix
H = alpha * Ham / sqrt(N);

%%% Step 4
P = zeros(N^2);
for k1=1:M
    P=P+D(k1,k1)/2*(2*kron(eye(N),DISS{k1})*kron(transpose(DISS{k1}'),eye(N))-kron(transpose(DISS{k1}'*DISS{k1}),eye(N))-kron(eye(N),DISS{k1}'*DISS{k1}));
end
HamPart = -sqrt(-1)*(kron(eye(N),H)-kron(transpose(H),eye(N)));
P = P + HamPart;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Lindbladian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
