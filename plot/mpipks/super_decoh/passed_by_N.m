clear all;
addpath('../../../source/matlab/lib')

Ns = linspace(10, 100, 91)';
p = 1e-10;
seeds = linspace(1, 100, 100)';

evals_lim = 1e-8;

passed_evals = zeros(size(Ns, 1) * size(seeds, 1), 1);
xs_global = zeros(size(Ns, 1) * size(seeds, 1), 1);
means = zeros(size(Ns, 1), 1);

for N_id = 1:size(Ns, 1)
    
    N = Ns(N_id);
    fprintf('N = %d\n', N);
    
    path = '/data/condmat/ivanchen/yusipov/os_lnd/super_decoh/eigen_dense';
    figures_path = '/home/ivanchen/yusipov/os_lnd/figures/super_decoh';
    
    N2 = N * N;
    
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
        evals = N * (real(evals) + 1) + 1i * N * imag(evals);
        
        curr_passed_evals = 0;
        for e_id = 1:size(evals, 1)
            if abs(imag(evals(e_id))) > evals_lim
                curr_passed_evals = curr_passed_evals + 1;
            end
        end
        
        id = (N_id - 1) * size(seeds, 1) + seed_id;
        passed_evals(id) = curr_passed_evals / N;
		xs_global(id) = N;
    end
    s_id = (N_id - 1) * size(seeds, 1) + 1;
    f_id = N_id * size(seeds, 1);
    means(N_id) = mean(passed_evals(s_id : f_id));
end

suffix = sprintf('N(var)_p(%0.10f)_numSeeds(%d)_logLim(%0.4f)', p, size(seeds, 1), log10(evals_lim));

min_y = min(passed_evals, [], 'all');
max_y = max(passed_evals, [], 'all');
fprintf('min_y = %0.16e\n', min_y);
fprintf('max_y = %0.16e\n', max_y);

pdf2d.x_num_bins = size(Ns, 1);
pdf2d.y_num_bins = 100;
pdf2d.x_label = '$N$';
pdf2d.y_label = '$passed/N$';
pdf2d.x_bin_s = Ns(1) - 0.5;
pdf2d.x_bin_f = Ns(end) + 0.5;
pdf2d.y_bin_s = min_y;
pdf2d.y_bin_f = max_y;
pdf2d = oqs_pdf_2d_setup(pdf2d);
data2d = horzcat(xs_global, passed_evals);
pdf2d = oqs_pdf_2d_update(pdf2d, data2d);
pdf2d = oqs_pdf_2d_release(pdf2d);
fig = oqs_pdf_2d_plot(pdf2d);
hold all;
plot(Ns, means, 'LineWidth', 3);
fn_fig = sprintf('%s/passed_by_N_%s', figures_path, suffix);
oqs_save_fig(fig, fn_fig)
