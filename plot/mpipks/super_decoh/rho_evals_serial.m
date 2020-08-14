clear all;
addpath('../../../source/matlab/lib')

figures_path = '/home/ivanchen/yusipov/os_lnd/figures/super_decoh';

task = 'eigen_dense';
path = sprintf('/data/condmat/ivanchen/yusipov/os_lnd/serial/super_decoh/%s', task);

method = 'simple';
G_type = 0;
ad = 0;
reshuffle_type = 0;

Ns = [64]';
ps = [1]';

num_seeds_total = 100000;

for N_id = 1:size(Ns, 1)
    
    N = Ns(N_id);
    N2 = N * N;
    
    for p_id = 1:size(ps, 1)
        
        p = ps(p_id);
        fprintf('p = %0.16e\n', p);
        
        rho_evals_all = zeros(N * num_seeds_total, 1);
        
        prefix = sprintf('method_%s/G_type_%d_ad_%d/reshuffle_type_%d/N_%d/p_%0.4f', ...
            method, ...
            G_type, ...
            ad, ...
            reshuffle_type, ...
            N, ...
            p);
        
        aggr_path = sprintf('%s/aggregator/%s', ...
            path, ...
            prefix);
        mkdir(aggr_path);
        
        suffix = sprintf('reshuffle(%d)_G(%d)_ad(0)_N(%d)_p(%0.4f)_seeds(%d_%d_%d)', ...
            reshuffle_type, ...
            G_type, ...
            N, ...
            p, ...
            1, ...
            1, ...
            num_seeds_total);
		
		fn = sprintf('%s/rho_evals_%s.txt', aggr_path, suffix)
		 
		rho_evals_data = importdata(fn);
		all_evals = rho_evals_data(:, 1);
		
		evals_data.x_num_bins = 400;
		evals_data.x_bin_s = min(all_evals);
		evals_data.x_bin_f = max(all_evals);
		evals_data.x_label = '$\lambda$';
		
		evals_data = oqs_pdf_1d_setup(evals_data);
		evals_data = oqs_pdf_1d_update(evals_data, all_evals);
		evals_data = oqs_pdf_1d_release(evals_data);
		fig = oqs_pdf_1d_plot(evals_data);
		fn_fig = sprintf('%s/rho_evals_%s', figures_path, suffix);
		oqs_save_fig(fig, fn_fig)
		
		mean_evals = mean(all_evals);
		fprintf('mean_evals = %0.16e\n', mean_evals);
		all_evals_moved = all_evals - mean_evals;
		
		evals_data_moved.x_num_bins = 400;
		evals_data_moved.x_bin_s = min(all_evals_moved);
		evals_data_moved.x_bin_f = max(all_evals_moved);
		evals_data_moved.x_label = '$\lambda$';
		evals_data_moved = oqs_pdf_1d_setup(evals_data_moved);
		evals_data_moved = oqs_pdf_1d_update(evals_data_moved, all_evals_moved);
		evals_data_moved = oqs_pdf_1d_release(evals_data_moved);
		fig = oqs_pdf_1d_plot(evals_data_moved);
		fn_fig = sprintf('%s/rho_evals_moved_%s', figures_path, suffix);
		oqs_save_fig(fig, fn_fig)
    end
end
