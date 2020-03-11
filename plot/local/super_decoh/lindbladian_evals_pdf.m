clear all;

N = 50;
p = 0.01;
seed = 10;

pdf2d.x_num_bins = 101;
pdf2d.y_num_bins = 101;
pdf2d.x_label = '$Re(\lambda)$';
pdf2d.y_label = '$Im(\lambda)$';

pdf1dlog.x_num_bins = 101;
pdf1dlog.x_label = '$Im(\lambda)$';

path = '../../../source/cpp/os_lnd/os_lnd';
suffix = sprintf('N(%d)_p(%0.4f)_seed(%d)', N, p, seed);

N2 = N * N;

evals = zeros(N2, 1);
fn_cpp = sprintf('%s/lindbladian_evals_%s.txt', path, suffix);
cpp_data = importdata(fn_cpp);
for str_id = 1:N2
    str = string(cpp_data(str_id));
    data2d = sscanf(str, '(%e,%e)', 2);
    evals(str_id) = data2d(1) + 1i * data2d(2);
end
evals = sort(evals, 'ComparisonMethod', 'abs');
evals = evals(2:end);
evals = N * (real(evals) + 1) + 1i * N * imag(evals);

abs_imag_parts = sort(abs(imag(evals)), 'descend');
pdf1dlog.x_bin_s = min(abs_imag_parts);
pdf1dlog.x_bin_f = max(abs_imag_parts);
pdf1dlog = oqs_pdf_1d_log_setup(pdf1dlog);
pdf1dlog = oqs_pdf_1d_log_update(pdf1dlog, abs_imag_parts);
pdf1dlog = oqs_pdf_1d_log_release(pdf1dlog);
fig = oqs_pdf_1d_log_plot(pdf1dlog);
fn_fig = sprintf('abs_imag_parts_%s', suffix);
oqs_save_fig(fig, fn_fig)

pdf1dlog.x_bin_s = min(abs_imag_parts);
pdf1dlog.x_bin_f = max(abs_imag_parts);



pdf2d.x_bin_s = min(real(evals));
pdf2d.x_bin_f = max(real(evals));
pdf2d.y_bin_s = min(imag(evals));
pdf2d.y_bin_f = max(imag(evals));
pdf2d = oqs_pdf_2d_setup(pdf2d);
data2d = horzcat(real(evals), imag(evals));
pdf2d = oqs_pdf_2d_update(pdf2d, data2d);
pdf2d = oqs_pdf_2d_release(pdf2d);
fig = oqs_pdf_2d_plot(pdf2d);
fn_fig = sprintf('lindbladian_evals_%s', suffix);
oqs_save_fig(fig, fn_fig)
