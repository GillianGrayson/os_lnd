clear all;

addpath('../../../source/matlab/lib')

task = 'eigen_dense';
path = sprintf('/data/biophys/denysov/yusipov/os_lnd/lind_ham/%s', task);
figures_path = '/home/denysov/yusipov/os_lnd/figures/lind_ham';

imag_lim = 5e-3;

alpha = 0.5;
N = 100; % system size
num_seeds = 100;

M = N^2-1; % auxiliary size

seeds = linspace(1, num_seeds, num_seeds)';


all_evals = zeros(M * num_seeds, 1);
zs_all = zeros(M * num_seeds, 1);
s_id = 0;
f_id = 0;
for seed = 1:num_seeds
    
    prefix = sprintf('N_%d/alpha_%0.4f/seed_%d', ...
        N, ...
        alpha, ...
        seed);
    
    suffix = sprintf('gen(1)_N(%d)_alpha(%0.4f)_seed(%d)', ...
        N, ...
        alpha, ...
        seed);
    
    fn = sprintf('%s/lindbladian_evals_%s.txt', aggr_path, suffix);

    evals = importdata(fn);
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
fn_fig = sprintf('%s/zs_alpha(%0.2f)_N(%d)_num_seeds(%d)', figures_path, alpha, N, num_seeds);
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
fn_fig = sprintf('%s/Levals_alpha(%0.2f)_N(%d)_num_seeds(%d)', figures_path, alpha, N, num_seeds);
oqs_save_fig(fig, fn_fig)

