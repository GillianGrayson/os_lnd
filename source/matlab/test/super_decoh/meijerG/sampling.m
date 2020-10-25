clear all;

N = 100;
num_seeds = 10;
N2 = N^2;

pdf.x_num_bins = 200;
pdf.x_label = '$\lambda$';
    
all_seeds = zeros(N2, num_seeds);

for seed = 1:num_seeds
    fprintf('seed = %d\n', seed);
    C = Hgen(seed, N);
    B = kron(eye(N), C) + kron(C, eye(N));
    evals = eig(B);
    all_seeds(:, seed) = evals;

end
pdf.x_bin_s = min(all_seeds, [], 'all')
pdf.x_bin_f = max(all_seeds, [], 'all')
pdf = oqs_pdf_1d_setup(pdf);

for seed = 1:num_seeds
    pdf = oqs_pdf_1d_update(pdf, all_seeds(:, seed));
end

pdf = oqs_pdf_1d_release(pdf);
fig = oqs_pdf_1d_plot(pdf);
fn_fig = sprintf('sampling');
oqs_save_fig(fig, fn_fig);
