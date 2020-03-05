N = 10; % system size
p = 0.1; % parameter for partial decoherence 0<=p<=1
seed = 10;

check_f_basis = 0;

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

FF1 = Reshuf(P,N,N2);

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


L_gamma = zeros(N2);
for i = 1:N
    for j = 1:N
        for k = 1:N
            for m = 1:N
                origin_row = (i - 1) * N + k;
                origin_col = (j - 1) * N + m;
                
                gamma_row = (i - 1) * N + j;
                gamma_col = (k - 1) * N + m;
                
                L_gamma(gamma_row, gamma_col) = P(origin_row, origin_col);
            end
        end
    end
end

P_gamma = zeros(N2);
for i = 1:N
    for j = 1:N
        for k = 1:N
            for m = 1:N
                origin_row = (i - 1) * N + k;
                origin_col = (j - 1) * N + m;
                
                gamma_row = (m - 1) * N + k;
                gamma_col = (j - 1) * N + i;
                
                P_gamma(gamma_row, gamma_col) = P(origin_row, origin_col);
            end
        end
    end
end

lnd_diff = norm(L_gamma - L)
lnd_diff = norm(P_gamma - FF1)

lnd_diff = norm(P_gamma - L_gamma)

 % transfroming the map into a state (non-normalized)
FF = Decoh(FF1,N2,p); % decoherence acting on the state
trace(FF1*FF1)
trace(FF1)
P = Reshuf(FF,N,N2); % transforming the state back into a map

% caclualting an addition to the map in order to get the Lindbladian
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

trace(AA)

%%%%% generation of Lindblad operator L %%%%%%
P = P - (1/2)*(kron(AA,eye(N)) + kron(eye(N),transpose(AA)));
%%%%%
% end of formaing Lindbladian

%%% END OF SUPERDECOHERENCE ACTION


%Step 5
[Evec,D] = eig(P); % eigenvalues of Lindblad superoperator
Eval = diag(D);

Eval_re = real(Eval);
Eval_im = imag(Eval);

fname1 = strcat('eval_L_re_',str1,'.txt');
fname2 = strcat('eval_L_im_',str1,'.txt');

hull = fopen('sobstvenie_znacheniya.dat','w');
%fprintf(hull,'%8.4f %16.12f\n',DT.Points(C,1),DT.Points(C,2));
for nn=2:M
    sx=N*real(Eval(nn)+1);
    sy=N*imag(Eval(nn));
    fprintf(hull,'%12.8f %16.12f\n',sx,sy);
end







%write evals of L-matrix to txt-files
save(fname1, 'Eval_re', '-ASCII', '-double')
save(fname2, 'Eval_im', '-ASCII', '-double')

%{
rho = zeros(N,N);

for i=1:N
    for j=1:N
        rho(j,i) = Evec((i-1)*N+j,1);
    end
end

rho = rho/trace(rho);

rho_re = real(rho);
rho_im = imag(rho);
%}

%{
fname3 = strcat('rho_re_',str1,'.txt');
fname4 = strcat('rho_im_',str1,'.txt');

%write rho-matrix to txt-files
save(fname3, 'rho_re', '-ASCII', '-double')
save(fname4, 'rho_im', '-ASCII', '-double')


rhoEval = eig(rho);
rhoEval_re = real(rhoEval);

fname5 = strcat('rhoEval_',str1,'.txt');
%write eigenvalues of rho-matrix to txt-files
save(fname5, 'rhoEval_re', '-ASCII', '-double')
%}

%quit()
