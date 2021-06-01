clear all;

N = 100;
num_seeds = 1000;
N2 = N^2;

pdf.x_num_bins = 300;
pdf.x_label = '$\lambda$';

pdf_H.x_num_bins = 200;
pdf_H.x_label = '$\lambda$';

H_evals = zeros(N, num_seeds);
all_evals = zeros(N2, num_seeds);

for seed = 1:num_seeds
    tic
    fprintf('seed = %d\n', seed);
    C = Hgen(seed, N);
    H_evals(:, seed) = eig(C);
    
    B = kron(eye(N), C) + kron(C, eye(N));
    evals = eig(B);
    all_evals(:, seed) = evals;
    toc
end

pdf_H.x_bin_s = min(H_evals, [], 'all');
pdf_H.x_bin_f = max(H_evals, [], 'all');
pdf_H = oqs_pdf_1d_setup(pdf_H);

for seed = 1:num_seeds
    pdf_H = oqs_pdf_1d_update(pdf_H, H_evals(:, seed));
end

pdf_H = oqs_pdf_1d_release(pdf_H);
fig = oqs_pdf_1d_plot(pdf_H);
fn_fig = sprintf('sampling_H');
oqs_save_fig(fig, fn_fig);

pdf.x_bin_s = min(all_evals, [], 'all');
pdf.x_bin_f = max(all_evals, [], 'all');
pdf = oqs_pdf_1d_setup(pdf);

for seed = 1:num_seeds
    pdf = oqs_pdf_1d_update(pdf, all_evals(:, seed));
end

pdf = oqs_pdf_1d_release(pdf);
fig = oqs_pdf_1d_plot(pdf);
fn_fig = sprintf('sampling_B');
oqs_save_fig(fig, fn_fig);
