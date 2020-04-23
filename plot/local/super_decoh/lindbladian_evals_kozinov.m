clear all;

reshufle_type = 1;
G_type = 0;
N = 200;
ps = [1 0.8 0.5 0.3 0.2 0.1 0.01]';
scaling_types = [1 1 1 1 1 1 1]';
seeds = linspace(1, 1, 1)';

path = 'E:/YandexDisk/Work/os_lnd/super_decoh/N_200_scale.txt';
cpp_data = importdata(path);

for p_id = 1:size(ps, 1)

	p = ps(p_id);

	fprintf('p = %0.16e\n', p);
    
    pdf2d.x_num_bins = 201;
    pdf2d.y_num_bins = 201;
    pdf2d.x_label = '$Re(\lambda)$';
    pdf2d.y_label = '$Im(\lambda)$';
    
    pdf2d_rem.x_num_bins = 201;
    pdf2d_rem.y_num_bins = 201;
    pdf2d_rem.x_label = '$Re(\lambda)$';
    pdf2d_rem.y_label = '$Im(\lambda)$';
    
    path = 'E:/YandexDisk/Work/os_lnd/super_decoh';
    
    N2 = N * N;
    
    all_evals = zeros((N2 - 1) * size(seeds, 1), 1);
    all_evals_rem = zeros((N - 1) * size(seeds, 1), 1);
    all_num_passed_evals = zeros(size(seeds, 1), 1);
    all_evec_sub_diag_norms = zeros((N - 1) * size(seeds, 1), 1);
        
    evals = cpp_data(:, 2 * p_id - 1) + 1i * cpp_data(:, 2 * p_id);
    
    [evals, order] = sort(evals, 'ComparisonMethod', 'abs');
    evals = evals(1:end);
%     if (scaling_types(p_id) == 1)
%         evals = N * sqrt(N) * (real(evals) + 1) + 1i * N * sqrt(N) * imag(evals);
%     elseif (scaling_types(p_id) == 3)
%         evals = N / p * (real(evals) + 1) + 1i * N / p * imag(evals);
%     else
%         evals = N * (real(evals) + 1) + 1i * N * imag(evals);
%     end
    
    s_id = 1;
    f_id = (N2 - 1);
    all_evals(s_id : f_id) = evals;

    suffix = sprintf('reshuffle(%d)_N(%d)_p(%0.10f)_numSeeds(%d)', reshufle_type, N, p, size(seeds, 1));
    
    pdf2d.x_bin_s = -4;
    pdf2d.x_bin_f = 4;
    pdf2d.y_bin_s = -1;
    pdf2d.y_bin_f = 1;
    pdf2d = oqs_pdf_2d_setup(pdf2d);
    data2d = horzcat(real(all_evals), imag(all_evals));
    pdf2d = oqs_pdf_2d_update(pdf2d, data2d);
    pdf2d = oqs_pdf_2d_release(pdf2d);
    fig = oqs_pdf_2d_plot(pdf2d);
    fn_fig = sprintf('lindbladian_evals_%s', suffix);
    oqs_save_fig(fig, fn_fig)
    
    pdf2d_rem.x_bin_s = -4;
    pdf2d_rem.x_bin_f = 4;
    pdf2d_rem.y_bin_s = -1;
    pdf2d_rem.y_bin_f = 1;
    pdf2d_rem = oqs_pdf_2d_setup(pdf2d_rem);
    data2d = horzcat(real(all_evals_rem), imag(all_evals_rem));
    pdf2d_rem = oqs_pdf_2d_update(pdf2d_rem, data2d);
    pdf2d_rem = oqs_pdf_2d_release(pdf2d_rem);
    fig = oqs_pdf_2d_plot(pdf2d_rem);
    fn_fig = sprintf('lindbladian_evals_rem_%s', suffix);
    oqs_save_fig(fig, fn_fig)
end