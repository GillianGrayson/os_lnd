clear all;

N = 10; % system size
p = 0.1; % parameter for partial decoherence 0<=p<=1
seed = 1;

check_f_basis = 1;
reshufle_type = 0; % 0 - original, 1 - mine

cpp_path = 'E:/Work/os_lnd/source/cpp/os_lnd/os_lnd';
suffix = sprintf('N(%d)_p(%0.4f)_seed(%d)', N, p, seed);

N2 = N^2;
M = N^2-1; % auxiliary size
imag1 = sqrt(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             f_basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=1;
F{1}=sparse(eye(N))/sqrt(N);
for i=1:N
    for j=i+1:N
        k=k+1;
        F{k}=sparse([i j],[j i],[1 1]/sqrt(2),N,N);
        k=k+1;
        F{k}=sparse([i j],[j i],-imag1*[1 -1]/sqrt(2),N,N);
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
%                            Lindbladian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=zeros(N^2);
for k1=1:M
    for k2=1:M
        %P=P+G(k1,k2)/2*(2*kron(eye(N),F{k1+1})*kron(transpose(F{k2+1}'),eye(N)));
        P=P+G(k1,k2)*kron(F{k1+1},eye(N))*kron(eye(N),transpose(F{k2+1}));
    end
end

if reshufle_type == 1
    FF1 = reshuffle(P, N);
    FF = Decoh(FF1, N2, p);
    P = reshuffle(FF, N);
else
    FF1 = Reshuf(P,N,N2);
    FF = Decoh(FF1,N2,p);
    P = Reshuf(FF,N,N2);
end

AA=zeros(N);

for s1 = 1:N
    for s2 = 1:N
        for s3 = 1:N
            w1=s3+N*(s1-1);
            w2=s3+N*(s2-1);
            AA(s1,s2)=AA(s1,s2)+ FF(w1,w2);
        end
    end
end

A = zeros(N);
fn_cpp = sprintf('%s/A_mtx_%s.txt', cpp_path, suffix);
cpp_data = importdata(fn_cpp);
for st_id_1 = 1:N
    for st_id_2 = 1:N
        str_id = (st_id_1 - 1) * N + st_id_2;
        str = string(cpp_data(str_id));
        data = sscanf(str, '(%e,%e)', 2);
        A(st_id_1, st_id_2) = data(1) + 1i * data(2);
    end
end
A = transpose(A); % Eigen is column-major by default.Here we have row-major dense matrix

a_diff = norm(AA - A)

P = P - (1/2)*(kron(AA,eye(N)) + kron(eye(N),transpose(AA)));

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               passed_evals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
passed_evals = [];
for evec_id = 1:N2
    evec_mtx = zeros(N,N);
    
    for i=1:N
        for j=1:N
            evec_mtx(i,j) = Evec((i-1)*N+j, evec_id);
        end
    end
    
    test_mtx = evec_mtx - diag(diag(evec_mtx));
    if norm(test_mtx) < 1e-12
        passed_evals = vertcat(passed_evals, Eval(evec_id));
    end

end

ololo = 1;


