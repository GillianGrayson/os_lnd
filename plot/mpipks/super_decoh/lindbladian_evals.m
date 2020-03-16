clear all;
addpath('../../../source/matlab/lib')

N = 100;
ps = logspace(-10, 0, 11);
seeds = linspace(1, 100, 100)';

evals_lim = 1e-8;

for p = ps

	fprintf('p = %0.16e\n', p);
    
    pdf2d.x_num_bins = 501;
    pdf2d.y_num_bins = 501;
    pdf2d.x_label = '$Re(\lambda)$';
    pdf2d.y_label = '$Im(\lambda)$';
    
    pdf2d_paased.x_num_bins = 101;
    pdf2d_paased.y_num_bins = 101;
    pdf2d_paased.x_label = '$Re(\lambda)$';
    pdf2d_paased.y_label = '$Im(\lambda)$';
    
    pdf1dlog.x_num_bins = 101;
    pdf1dlog.x_label = '$Im(\lambda)$';
    
    passed.x_num_bins = 5;
    passed.x_label = 'passed evals';
    
    path = '/data/condmat/ivanchen/yusipov/os_lnd/super_decoh/eigen_dense';
    figures_path = '/home/ivanchen/yusipov/os_lnd/figures/super_decoh';
    
    N2 = N * N;
    
    all_evals = zeros((N2 - 1) * size(seeds, 1), 1);
    all_evals_passed = [];
    num_passed_evals = zeros(size(seeds, 1), 1);
    
    for seed_id = 1:size(seeds, 1)
        seed = seeds(seed_id);
        
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
        evals = N * sqrt(N) * (real(evals) + 1) + 1i * N * sqrt(N) * imag(evals);
        
        curr_passed_evals = 0;
        for e_id = 1:size(evals, 1)
            if abs(imag(evals(e_id))) > evals_lim
                curr_passed_evals = curr_passed_evals + 1;
                all_evals_passed = vertcat(all_evals_passed, evals(e_id));
            end
        end
        num_passed_evals(seed_id) = curr_passed_evals;
        
        curr_passed_evals = curr_passed_evals;
        
        s_id = (seed_id - 1) * (N2 - 1) + 1;
        f_id = seed_id * (N2 - 1);
        all_evals(s_id : f_id) = evals;
    end
    
    suffix = sprintf('N(%d)_p(%0.10f)_numSeeds(%d)_logLim(%0.4f)', N, p, size(seeds, 1), log10(evals_lim));
    
    passed.x_bin_s = min(num_passed_evals);
    passed.x_bin_f = max(num_passed_evals);
    passed = oqs_pdf_1d_setup(passed);
    passed = oqs_pdf_1d_update(passed, num_passed_evals);
    passed = oqs_pdf_1d_release(passed);
    fig = oqs_pdf_1d_plot(passed);
    fn_fig = sprintf('%s/passed_evals_%s', figures_path, suffix);
    oqs_save_fig(fig, fn_fig);
    
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
    
    pdf2d_paased.x_bin_s = min(real(all_evals_passed));
    pdf2d_paased.x_bin_f = max(real(all_evals_passed));
    pdf2d_paased.y_bin_s = min(imag(all_evals_passed));
    pdf2d_paased.y_bin_f = max(imag(all_evals_passed));
    pdf2d_paased = oqs_pdf_2d_setup(pdf2d_paased);
    data2d = horzcat(real(all_evals_passed), imag(all_evals_passed));
    pdf2d_paased = oqs_pdf_2d_update(pdf2d_paased, data2d);
    pdf2d_paased = oqs_pdf_2d_release(pdf2d_paased);
    fig = oqs_pdf_2d_plot(pdf2d_paased);
    fn_fig = sprintf('%s/lindbladian_evals_passed_%s', figures_path, suffix);
    oqs_save_fig(fig, fn_fig)
    
end