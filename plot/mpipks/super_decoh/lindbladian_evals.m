clear all;
addpath('../../../source/matlab/lib')

method = 'origin';
reshufle_type = 1;
G_type = 1;
aux_dim = 5000;
N = 100;
ps = [1]';
scaling_types = [0]';
seeds = linspace(1, 100, 100)';

evals_lim = 1e-8;

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
    
    pdf1dlog.x_num_bins = 201;
    pdf1dlog.x_label = '$Im(\lambda)$';
    
    num_passed_evals.x_num_bins = 5;
    num_passed_evals.x_label = 'passed evals';

    evec_sub_diag_norms.x_num_bins = 50;
    evec_sub_diag_norms.x_label = '$\log_{10}norm$';
    
    path = '/data/condmat/ivanchen/yusipov/os_lnd/super_decoh/eigen_dense';
    figures_path = '/home/ivanchen/yusipov/os_lnd/figures/super_decoh';
    
    N2 = N * N;
    
    all_evals = zeros((N2 - 1) * size(seeds, 1), 1);
    all_evals_rem = zeros((N - 1) * size(seeds, 1), 1);
    all_num_passed_evals = zeros(size(seeds, 1), 1);
    all_evec_sub_diag_norms = zeros((N - 1) * size(seeds, 1), 1);
    
    for seed_id = 1:size(seeds, 1)
        
        seed = seeds(seed_id)
        
		suffix = sprintf('reshuffle(%d)_G(%d)_N(%d)_ad(%d)_p(%0.10f)_seed(%d)', reshufle_type, G_type, N, aux_dim, p, seed);
		
        evals = zeros(N2, 1);
		fn_cpp = sprintf('%s/method_%s/G_type_%d_ad_%d/reshuffle_type_%d/N_%d/p_%0.10f/seed_%d/lindbladian_evals_%s.txt', path, method, G_type, aux_dim, reshufle_type, N, p, seed, suffix);
        if ~isfile(fn_cpp)
			fprintf('Warning! seed = %d\n', seed);
		end
		cpp_data = importdata(fn_cpp);
        for str_id = 1:N2
            str = string(cpp_data(str_id));
            data2d = sscanf(str, '(%e,%e)', 2);
            evals(str_id) = data2d(1) + 1i * data2d(2);
        end
		
        evecs_norms = zeros(N2, 1);
        fn_cpp = sprintf('%s/method_%s/G_type_%d_ad_%d/reshuffle_type_%d/N_%d/p_%0.10f/seed_%d/evec_sub_diag_norms_%s.txt', path, method, G_type, aux_dim, reshufle_type, N, p, seed, suffix);
        cpp_data = importdata(fn_cpp);
        for str_id = 1:N2
            str = string(cpp_data(str_id));
            data1d = sscanf(str, '%e', 1 );
            evecs_norms(str_id) = data1d(1);
        end
        
        [evals, order] = sort(evals, 'ComparisonMethod', 'abs');
        evecs_norms = evecs_norms(order);
        evals = evals(2:end);
		if (scaling_types(p_id) == 1)
			evals = N * sqrt(N) * (real(evals) + 1) + 1i * N * sqrt(N) * imag(evals);
        elseif (scaling_types(p_id) == 3)
			evals = N / p * (real(evals) + 1) + 1i * N / p * imag(evals);
        elseif (scaling_types(p_id) == 0)
            nothing_to_do_here = 0;
		else
			evals = N * (real(evals) + 1) + 1i * N * imag(evals);
		end
		evecs_norms = evecs_norms(2:end);
        
        [evecs_norms, order] = sort(evecs_norms);
        
        for e_id = 1:N-1
            curr_norm = log10(evecs_norms(e_id) + 1e-16);
            all_evec_sub_diag_norms((seed_id - 1) * (N - 1) + e_id) = curr_norm;
            curr_eval = evals(order(e_id));
            all_evals_rem((seed_id - 1) * (N - 1) + e_id) = curr_eval;
        end
        
        curr_num_passed_evals = 0;
        for e_id = 1:size(evals, 1)
            if abs(imag(evals(e_id))) > evals_lim
                curr_num_passed_evals = curr_num_passed_evals + 1;
            end
        end
        all_num_passed_evals(seed_id) = curr_num_passed_evals;
        
        s_id = (seed_id - 1) * (N2 - 1) + 1;
        f_id = seed_id * (N2 - 1);
        all_evals(s_id : f_id) = evals;
    end
    
    suffix = sprintf('method(%s)_reshuffle(%d)_N(%d)_ad(%d)_p(%0.10f)_numSeeds(%d)_logLim(%0.4f)', method, reshufle_type, N, aux_dim, p, size(seeds, 1), log10(evals_lim));
    
    num_passed_evals.x_bin_s = min(all_num_passed_evals);
    num_passed_evals.x_bin_f = max(all_num_passed_evals);
    num_passed_evals = oqs_pdf_1d_setup(num_passed_evals);
    num_passed_evals = oqs_pdf_1d_update(num_passed_evals, all_num_passed_evals);
    num_passed_evals = oqs_pdf_1d_release(num_passed_evals);
    fig = oqs_pdf_1d_plot(num_passed_evals);
    fn_fig = sprintf('%s/passed_evals_%s', figures_path, suffix);
    oqs_save_fig(fig, fn_fig);
    
    suffix = sprintf('method(%s)_reshuffle(%d)_N(%d)_ad(%d)_p(%0.10f)_numSeeds(%d)', method, reshufle_type, N, aux_dim, p, size(seeds, 1));
    
    abs_imag_parts = sort(abs(imag(all_evals)), 'descend');
    pdf1dlog.x_bin_s = max(min(abs_imag_parts), 1e-17);
    pdf1dlog.x_bin_f = max(abs_imag_parts);
    pdf1dlog = oqs_pdf_1d_log_setup(pdf1dlog);
    pdf1dlog = oqs_pdf_1d_log_update(pdf1dlog, abs_imag_parts);
    pdf1dlog = oqs_pdf_1d_log_release(pdf1dlog);
    fig = oqs_pdf_1d_log_plot(pdf1dlog);
    fn_fig = sprintf('%s/abs_imag_parts_%s', figures_path, suffix);
    oqs_save_fig(fig, fn_fig);
    
    if (scaling_types(p_id) == 1)
        pdf2d.x_bin_s = -4;
        pdf2d.x_bin_f = 4;
        pdf2d.y_bin_s = -1;
        pdf2d.y_bin_f = 1;
    elseif (scaling_types(p_id) == 3)
        pdf2d.x_bin_s = -2;
        pdf2d.x_bin_f = 2;
        pdf2d.y_bin_s = -1;
        pdf2d.y_bin_s = 1;
    elseif (scaling_types(p_id) == 0)
        pdf2d.x_bin_s = min(real(all_evals));
        pdf2d.x_bin_f = max(real(all_evals));
        pdf2d.y_bin_s = min(imag(all_evals));
        pdf2d.y_bin_f = max(imag(all_evals));
    end
    pdf2d = oqs_pdf_2d_setup(pdf2d);
    data2d = horzcat(real(all_evals), imag(all_evals));
    pdf2d = oqs_pdf_2d_update(pdf2d, data2d);
    pdf2d = oqs_pdf_2d_release(pdf2d);
    fig = oqs_pdf_2d_plot(pdf2d);
    
    limit_pdf = 0.00000001;
    xr1 = pdf2d.x_bin_centers';
    yr1 = zeros(size(pdf2d.x_bin_centers, 2), 1);
    for x_id = 1:size(pdf2d.x_bin_centers, 2)
        is_found = 0;
        for y_id = 1:size(pdf2d.y_bin_centers, 2)
            if pdf2d.pdf(x_id, y_id) > limit_pdf
                yr1(x_id) = pdf2d.y_bin_centers(y_id);
                is_found = 1;
                break;
            end
            if is_found == 0
                yr1(x_id) = 0;
            end
        end
    end
    yr1 = smooth(yr1 * -1.0);
    xr2 = flip(xr1);
    yr2 = -flip(yr1);
    xr = vertcat(xr1, xr2);
    yr = vertcat(yr1, yr2);
    
    figure(fig)
    hold all;
    plot(xr, yr, 'LineWidth', 2)
    
    fn_fig = sprintf('%s/lindbladian_evals_%s', figures_path, suffix);
    oqs_save_fig(fig, fn_fig)
    
    pdf2d_rem.x_bin_s = pdf2d.x_bin_s;
    pdf2d_rem.x_bin_f = pdf2d.x_bin_f;
    pdf2d_rem.y_bin_s = pdf2d.y_bin_s;
    pdf2d_rem.y_bin_f = pdf2d.y_bin_s;
    pdf2d_rem = oqs_pdf_2d_setup(pdf2d_rem);
    data2d = horzcat(real(all_evals_rem), imag(all_evals_rem));
    pdf2d_rem = oqs_pdf_2d_update(pdf2d_rem, data2d);
    pdf2d_rem = oqs_pdf_2d_release(pdf2d_rem);
    fig = oqs_pdf_2d_plot(pdf2d_rem);
    fn_fig = sprintf('%s/lindbladian_evals_rem_%s', figures_path, suffix);
    oqs_save_fig(fig, fn_fig)
    
    evec_sub_diag_norms.x_bin_s = min(all_evec_sub_diag_norms) - 1e-16;
    evec_sub_diag_norms.x_bin_f = max(all_evec_sub_diag_norms) + 1e-16;
    evec_sub_diag_norms = oqs_pdf_1d_setup(evec_sub_diag_norms);
    evec_sub_diag_norms = oqs_pdf_1d_update(evec_sub_diag_norms, all_evec_sub_diag_norms);
    evec_sub_diag_norms = oqs_pdf_1d_release(evec_sub_diag_norms);
    fig = oqs_pdf_1d_plot(evec_sub_diag_norms);
    fn_fig = sprintf('%s/evec_sub_diag_norms_%s', figures_path, suffix);
    oqs_save_fig(fig, fn_fig);
    
end