clear all;

N = 5;
tau = 1;
k = 1;

con = zeros(4);
con(1, 1) = tau;
con(4, 4) = k;
con(2, 3) = 1;
con(3, 2) = 1;

dissipators = cell(N, 1);

left_diss = con;
for spin_id = 3:N
    left_diss = kron(left_diss, eye(2));
end
dissipators{1} = left_diss;

right_diss = eye(2);
for spin_id = 2:N-2
    right_diss = kron(right_diss, eye(2));
end
right_diss = kron(right_diss, con);
dissipators{N-1} = left_diss;


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

border_diss = zeros(2.^N);
num_states = size(border_diss, 1);
for state_id = 1:num_states
    state_int = state_id - 1;
    state_bits = de2bi(state_int, N);
    
    tmp = bi2de(state_bits);
    
    ololo = 1;
end
    

imag_lim = 5e-3;

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
    
all_evals = zeros(M * num_seeds, 1);
zs_all = zeros(M * num_seeds, 1);
s_id = 0;
f_id = 0;
for seed = 1:num_seeds
    fprintf('seed = %d\n', seed);
    tic
    
    %%% Step 2: calculate G-matrix
    G = Ggen(seed, M, N);
    
    [Evec,D] = eig(G);
    for k1=1:M
        DISS{k1}=0;
        for k2=1:M
            DISS{k1}=DISS{k1} + Evec(k2,k1) * F{k2+1};
        end
    end

    %%% Step 3: calculate H-matrix
    H = alpha * Hgen(seed, N) / sqrt(N);

    %%% Step 4
    P_origin = zeros(N^2); %% Initialize superopetator matrix
    P_diss = zeros(N^2);
    for k1=1:M
        
        P_diss=P_diss+D(k1,k1)/2*(2*kron(eye(N),DISS{k1})*kron(transpose(DISS{k1}'),eye(N))-kron(transpose(DISS{k1}'*DISS{k1}),eye(N))-kron(eye(N),DISS{k1}'*DISS{k1}));
        
        for k2=1:M
            P_origin=P_origin+G(k1,k2)/2*(2*kron(eye(N),F{k1+1})*kron(transpose(F{k2+1}'),eye(N))-kron(transpose(F{k2+1}'*F{k1+1}),eye(N))-kron(eye(N),F{k2+1}'*F{k1+1}));
        end
    end
    HamPart = -sqrt(-1)*(kron(eye(N),H)-kron(transpose(H),eye(N)));
    P_origin = P_origin + HamPart;
    P_diss = P_diss + HamPart;
    
    P_check = norm(P_origin - P_diss);
    if (abs(P_check) > 1e-12)
        fprintf("Lindbladians are differ!")
    end

    %%% Step 5
    evals = eig(P_origin);
    evals = sort(evals,'ComparisonMethod','abs');
    evals = evals(2:end);
    evals = (evals + 1) * N;
    all_evals((seed - 1) * M + 1 : seed * M) = evals;
    
    evals_filtered = evals(abs(imag(evals)) >= imag_lim);
    num_passed_evals = size(evals_filtered, 1);
    fprintf('num_passed_evals = %d\n', num_passed_evals);
    
    s_id = f_id + 1;
    f_id = s_id + num_passed_evals - 1;
    
    zs = zeros(num_passed_evals, 1);
    neighbours_2D = horzcat(real(evals_filtered), imag(evals_filtered));
    for z_id = 1 : num_passed_evals
        target_2D = horzcat(real(evals_filtered(z_id)), imag(evals_filtered(z_id)));
        target = evals_filtered(z_id);
        distances = sqrt((neighbours_2D(:, 1) - target_2D(:, 1)).^2 + (neighbours_2D(:, 2) - target_2D(:, 2)).^2);
        [min_distances, order] = mink(distances, 3);
        nn = evals_filtered(order(2));
        nnn = evals_filtered(order(3));
        zs(z_id) = (nn - target) / (nnn - target);
    end
    
    num_zero_zs = sum(abs(zs) < 1e-4);
    fprintf('num_zero_zs = %d\n', num_zero_zs);
    
    zs_all(s_id : f_id) = zs;
    
    toc
end

zs_all = zs_all(1 : f_id);

total_num_passed_evals = f_id

pdf2dzs.x_num_bins = 101;
pdf2dzs.y_num_bins = 101;
pdf2dzs.x_label = '$Re(z)$';
pdf2dzs.y_label = '$Im(z)$';

pdf2dzs.x_bin_s = -1;
pdf2dzs.x_bin_f = 1;
pdf2dzs.y_bin_s = -1;
pdf2dzs.y_bin_f = 1;
pdf2dzs = oqs_pdf_2d_setup(pdf2dzs);
data2d = horzcat(real(zs_all), imag(zs_all));
pdf2dzs = oqs_pdf_2d_update(pdf2dzs, data2d);
pdf2dzs = oqs_pdf_2d_release(pdf2dzs);
fig = oqs_pdf_2d_plot(pdf2dzs);
fn_fig = sprintf('zs_alpha(%0.2f)_N(%d)_num_seeds(%d)', alpha, N, num_seeds);
oqs_save_fig(fig, fn_fig)

pdf2d.x_num_bins = 201;
pdf2d.y_num_bins = 201;
pdf2d.x_label = '$Re(\lambda)$';
pdf2d.y_label = '$Im(\lambda)$';

pdf2d.x_bin_s = min(real(all_evals));
pdf2d.x_bin_f = max(real(all_evals));
pdf2d.y_bin_s = min(imag(all_evals));
pdf2d.y_bin_f = max(imag(all_evals));

pdf2d = oqs_pdf_2d_setup(pdf2d);
data2d = horzcat(real(all_evals), imag(all_evals));
pdf2d = oqs_pdf_2d_update(pdf2d, data2d);
pdf2d = oqs_pdf_2d_release(pdf2d);
fig = oqs_pdf_2d_plot(pdf2d);
fn_fig = sprintf('Levals_alpha(%0.2f)_N(%d)_num_seeds(%d)', alpha, N, num_seeds);
oqs_save_fig(fig, fn_fig)

