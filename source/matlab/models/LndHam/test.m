clear all;

alpha = 1.0;
N = 10; % system size
M = N^2-1; % auxiliary size
num_seeds = 500;

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

    %%% Step 3: calculate H-matrix
    H = alpha * Hgen(seed, N) / sqrt(N);

    %%% Step 4
    P=zeros(N^2); %% Initialize superopetator matrix
    for k1=1:M
        for k2=1:M
            P=P+G(k1,k2)/2*(2*kron(eye(N),F{k1+1})*kron(transpose(F{k2+1}'),eye(N))-kron(transpose(F{k2+1}'*F{k1+1}),eye(N))-kron(eye(N),F{k2+1}'*F{k1+1}));
        end
    end
    A = -sqrt(-1)*(kron(eye(N),H)-kron(transpose(H),eye(N)));
    P = P + A;

    %%% Step 5
    evals = eig(P);
    evals = sort(evals,'ComparisonMethod','abs');
    evals = evals(2:end);
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

