clear all;

seed = 1;
N = 10;

path = 'E:/Work/os_lnd/source/cpp/os_lnd/os_lnd';
suffix = sprintf('N(%d)_seed(%d)', N, seed);

M = N^2;
imag1 = sqrt(-1.);

MS = zeros(N);
fn_cpp = sprintf('%s/MS_mtx_%s.txt', path, suffix);
cpp_data = importdata(fn_cpp);
for st_id_1 = 1:N
    for st_id_2 = 1:N
        str_id = (st_id_1 - 1) * N + st_id_2;
        str = string(cpp_data(str_id));
        data = sscanf(str, '(%e,%e)', 2);
        MS(st_id_1, st_id_2) = data(1) + 1i * data(2);
    end
end
MS = transpose(MS); % Eigen is column-major by default. Here we have row-major dense matrix

MZ = zeros(N,N);
for j = 1:N
    for k = 1:N
        MZ(k,j) = MS(k,j)*conj(MS(k,j));
    end
end
MZ = sqrt(N)*MZ/trace(MZ);


X = zeros(M);
fn_cpp = sprintf('%s/X_mtx_%s.txt', path, suffix);
cpp_data = importdata(fn_cpp);
for st_id_1 = 1:M
    for st_id_2 = 1:M
        str_id = (st_id_1 - 1) * M + st_id_2;
        str = string(cpp_data(str_id));
        data = sscanf(str, '(%e,%e)', 2);
        X(st_id_1, st_id_2) = data(1) + 1i * data(2);
    end
end
X = transpose(X); % Eigen is column-major by default. Here we have row-major dense matrix

rho_temp=X*X';
for j = 1:N
    for k = 1:N
        s=N*(j-1)+k;
        
        ttt=rho_temp(s,s)/MZ(j,k);
        
        for sp = 1:M
            X(sp,s)=X(sp,s)/sqrt(ttt);
        end
        
    end
end


rho_temp=X*X';
rho_temp = sqrt(M)*rho_temp/trace(rho_temp);


rho=rho_temp;
%rho = Decoh(rho_temp,M,p); % decoherence of the random state
F = Reshuff(rho,N,M); % reshuffling procedure

AA=zeros(N);

for s1 = 1:N
    for s2 = 1:N
        for s3 = 1:N
            w1=s3+N*(s1-1);
            w2=s3+N*(s2-1);
            AA(s1,s2)=AA(s1,s2)+rho(w1,w2);
        end
    end
end

L = F - (1/2)*(kron(AA,eye(N)) + kron(eye(N),transpose(AA)));

L_cpp = zeros(M);
fn_cpp = sprintf('%s/lindbladian_mtx_%s.txt', path, suffix);
cpp_data = importdata(fn_cpp);
for st_id_1 = 1:M
    for st_id_2 = 1:M
        str_id = (st_id_1 - 1) * M + st_id_2;
        str = string(cpp_data(str_id));
        data = sscanf(str, '(%e,%e)', 2);
        L_cpp(st_id_1, st_id_2) = data(1) + 1i * data(2);
    end
end
L_cpp = transpose(L_cpp); % Eigen is column-major by default. Here we have row-major dense matrix

norm_diff_L = norm(L - L_cpp)

evals = eig(L); % calculation of Lindblad eigenvalues

