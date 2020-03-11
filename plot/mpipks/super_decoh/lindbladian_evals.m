clear all;
addpath('../../../source/matlab/lib')

N = 50;
p = 0.50;
seeds = linspace(1, 100, 100)';

pdf2d.x_num_bins = 501;
pdf2d.y_num_bins = 501;
pdf2d.x_label = '$Re(\lambda)$';
pdf2d.y_label = '$Im(\lambda)$';

pdf1dlog.x_num_bins = 101;
pdf1dlog.x_label = '$Im(\lambda)$';

path = '/data/condmat/ivanchen/yusipov/os_lnd/super_decoh/eigen_dense';
figures_path = '/home/ivanchen/yusipov/os_lnd/figures/super_decoh';

N2 = N * N;

all_evals = zeros((N2 - 1) * size(seeds, 1), 1);

for seed_id = 1:size(seeds, 1)
    seed = seeds(seed_id)
	
	suffix = sprintf('N(%d)_p(%0.10f)_seed(%d)', N, p, seed);
    
    evals = zeros(N2, 1);
    fn_cpp = sprintf('%s/N_%d/p_%0.10f/seed_%d/lindbladian_evals_%s.txt', path, N, p, seed, suffix);
    cpp_data = importdata(fn_cpp);
    for str_id = 1:N2
        str = string(cpp_data(str_id));
        data2d = sscanf(str, '(%e,%e)', 2);
        evals(str_id) = data2d(1) + 1i * data2d(2);
    end
    evals = sort(evals, 'ComparisonMethod', 'abs');
    evals = evals(2:end);
    evals = N * (real(evals) + 1) + 1i * N * imag(evals);
    
	s_id = (seed_id - 1) * (N2 - 1) + 1;
	f_id = seed_id * (N2 - 1);
    all_evals(s_id : f_id) = evals;
end

suffix = sprintf('N(%d)_p(%0.10f)_numSeeds(%d)', N, p, size(seeds, 1));

abs_imag_parts = sort(abs(imag(all_evals)), 'descend');
pdf1dlog.x_bin_s = max(min(abs_imag_parts), 1e-17);
pdf1dlog.x_bin_f = max(abs_imag_parts);
pdf1dlog = oqs_pdf_1d_log_setup(pdf1dlog);
pdf1dlog = oqs_pdf_1d_log_update(pdf1dlog, abs_imag_parts);
pdf1dlog = oqs_pdf_1d_log_release(pdf1dlog);
fig = oqs_pdf_1d_log_plot(pdf1dlog);
fn_fig = sprintf('%s/abs_imag_parts_%s', figures_path, suffix);
oqs_save_fig(fig, fn_fig);

pdf2d.x_bin_s = min(real(all_evals));
pdf2d.x_bin_f = max(real(all_evals));
pdf2d.y_bin_s = min(imag(all_evals));
pdf2d.y_bin_f = max(imag(all_evals));
pdf2d = oqs_pdf_2d_setup(pdf2d);
data2d = horzcat(real(all_evals), imag(all_evals));
pdf2d = oqs_pdf_2d_update(pdf2d, data2d);
pdf2d = oqs_pdf_2d_release(pdf2d);
fig = oqs_pdf_2d_plot(pdf2d);
fn_fig = sprintf('%s/lindbladian_evals_%s', figures_path, suffix);
oqs_save_fig(fig, fn_fig)
