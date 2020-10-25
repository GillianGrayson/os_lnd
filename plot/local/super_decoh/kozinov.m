clear all;

path_figures = 'E:/YandexDisk/Work/os_lnd/figures/super_decoh/N300_kozinov';
path_data = 'E:/YandexDisk/Work/os_lnd/data/N300_kozinov';
path_border_line = 'E:/YandexDisk/Work/os_lnd/data';

num_hits = 0

N = 300;
ps = [1]';
seeds = [linspace(0, 0, 1)]';

N2 = N * N;
num_seeds = size(seeds, 1);

theory_data_classical = importdata(sprintf('%s/borderline_classical.dat', path_border_line));
theory_data_quantum = importdata(sprintf('%s/borderline_quantum.dat', path_border_line));

iou_classical = zeros(size(ps, 1), 1);
iou_quantum = zeros(size(ps, 1), 1);
    
for p_id = 1:size(ps, 1)

	p = ps(p_id);

	fprintf('p = %0.16e\n', p);
    
    pdf2d_classical.x_num_bins = 201;
    pdf2d_classical.y_num_bins = 201;
    pdf2d_classical.x_label = '$Re(\lambda)$';
    pdf2d_classical.y_label = '$Im(\lambda)$';
    
    pdf2d_quantum.x_num_bins = 201;
    pdf2d_quantum.y_num_bins = 201;
    pdf2d_quantum.x_label = '$Re(\lambda)$';
    pdf2d_quantum.y_label = '$Im(\lambda)$';
    
    all_evals_classical = zeros((N2 - 1) * num_seeds, 1);
    all_evals_quantum = zeros((N2 - 1) * num_seeds, 1);
    
    for seed_id = 1:size(seeds, 1)
        seed = seeds(seed_id);
        suffix = sprintf('%0.6f__N(%d)_dt(0)_dp(0.0000)_g(1.0000)_T(1.0000)_seed(%d)', ...
            p, ...
            N, ...
            seed);
        fn = sprintf('%s/P_new_evals_%s.txt', path_data, suffix);
        evals_data = importdata(fn);
        
        evals = evals_data((seed_id - 1) * N2 + 1 : seed_id * (N2 - 1), 1)  + 1i * evals_data((seed_id - 1) * N2 + 1 : seed_id * (N2 - 1), 2);
        
        %[evals, ~] = sort(evals, 'ComparisonMethod', 'abs');     
        %evals = evals(2:end);
        evals_classical = N * sqrt(N) * (real(evals) + 1) + 1i * N * sqrt(N) * imag(evals);
        evals_quantum = N / p * (real(evals) + 1) + 1i * N / p * imag(evals);
        
        s_id = (seed_id - 1) * (N2 - 1) + 1;
        f_id = seed_id * (N2 - 1);
        all_evals_classical(s_id : f_id) = evals_classical;
        all_evals_quantum(s_id : f_id) = evals_quantum;
    end
    
    x_lim = max(abs(min(real(all_evals_classical))), abs(max(real(all_evals_classical))));
    pdf2d_classical.x_bin_s = -x_lim - 1e-16;
    pdf2d_classical.x_bin_f = x_lim + 1e-16;
    pdf2d_classical.y_bin_s = min(min(imag(all_evals_classical)) - 1e-16, -0.8);
    pdf2d_classical.y_bin_f = max(max(imag(all_evals_classical)) + 1e-16, 0.8);
    pdf2d_classical = oqs_pdf_2d_setup(pdf2d_classical);
    data2d = horzcat(real(all_evals_classical), imag(all_evals_classical));
    pdf2d_classical = oqs_pdf_2d_update(pdf2d_classical, data2d);
    fprintf('num_evals = %d\n', pdf2d_classical.inc_count);
    pdf2d_classical = oqs_pdf_2d_release(pdf2d_classical);
    fig_classical_1 = oqs_pdf_2d_plot(pdf2d_classical);
    fig_classical_2 = oqs_pdf_2d_plot(pdf2d_classical);
    hold all;
    xt1 = theory_data_classical(:, 1);
    yt1 = theory_data_classical(:, 2);
    xt2 = flip(xt1);
    yt2 = -flip(yt1);
    xt = vertcat(xt1, xt2);
    yt = vertcat(yt1, yt2);
    pgont_classical = polyshape(xt, yt);
    figure(fig_classical_1)
    hold all;
    hLine = plot(xt, yt, 'LineWidth', 3, 'Color', 'green');
    legend(hLine, 'classical');
    fn_fig_classical = sprintf('%s/lindbladian_evals_classical_N(%d)_p(%0.4f)', path_figures, N, p);
    oqs_save_fig(fig_classical_1, fn_fig_classical);
    
    x_lim = max(abs(min(real(all_evals_quantum))), abs(max(real(all_evals_quantum))));
    pdf2d_quantum.x_bin_s = min(-x_lim - 1e-16, -2);
    pdf2d_quantum.x_bin_f = max(x_lim + 1e-16, 2);
    pdf2d_quantum.y_bin_s = min(min(imag(all_evals_quantum)) - 1e-16, -1.1);
    pdf2d_quantum.y_bin_f = max(max(imag(all_evals_quantum)) + 1e-16, 1.1);
    pdf2d_quantum = oqs_pdf_2d_setup(pdf2d_quantum);
    data2d = horzcat(real(all_evals_quantum), imag(all_evals_quantum));
    pdf2d_quantum = oqs_pdf_2d_update(pdf2d_quantum, data2d);
    fprintf('num_evals = %d\n', pdf2d_quantum.inc_count);
    pdf2d_quantum = oqs_pdf_2d_release(pdf2d_quantum);
    fig_quantum_1 = oqs_pdf_2d_plot(pdf2d_quantum);
    fig_quantum_2 = oqs_pdf_2d_plot(pdf2d_quantum);
    hold all;
    xt1 = theory_data_quantum(:, 1);
    yt1 = theory_data_quantum(:, 2);
    xt2 = flip(xt1);
    yt2 = -flip(yt1);
    xt = vertcat(xt1, xt2);
    yt = vertcat(yt1, yt2);
    pgont_quantum = polyshape(xt, yt);
    figure(fig_quantum_1)
    hold all;
    hLine = plot(xt, yt, 'LineWidth', 3, 'Color', 'blue');
    legend(hLine, 'quantum');
    fn_fig_quantum = sprintf('%s/lindbladian_evals_quantum_N(%d)_p(%0.4f)', path_figures, N, p);
    oqs_save_fig(fig_quantum_1, fn_fig_quantum);
        
    limit_pdf_classical = num_hits / (size(all_evals_classical, 1) * pdf2d_classical.x_bin_shift * pdf2d_classical.y_bin_shift);
    limit_pdf_quantum = num_hits / (size(all_evals_quantum, 1) * pdf2d_quantum.x_bin_shift * pdf2d_quantum.y_bin_shift);
    
    xr1 = pdf2d_classical.x_bin_centers';
    yr1 = zeros(size(pdf2d_classical.x_bin_centers, 2), 1);
    for x_id = 1:size(pdf2d_classical.x_bin_centers, 2)
        is_found = 0;
        for y_id = 1:size(pdf2d_classical.y_bin_centers, 2)
            if pdf2d_classical.pdf(x_id, y_id) > limit_pdf_classical
                yr1(x_id) = pdf2d_classical.y_bin_centers(y_id);
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
    pgonr_classical = polyshape(xr, yr);
    
    xr1 = pdf2d_quantum.x_bin_centers';
    yr1 = zeros(size(pdf2d_quantum.x_bin_centers, 2), 1);
    for x_id = 1:size(pdf2d_quantum.x_bin_centers, 2)
        is_found = 0;
        for y_id = 1:size(pdf2d_quantum.y_bin_centers, 2)
            if pdf2d_quantum.pdf(x_id, y_id) > limit_pdf_quantum
                yr1(x_id) = pdf2d_quantum.y_bin_centers(y_id);
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
    pgonr_quantum = polyshape(xr, yr);
    
    poly_intersect_quantum = intersect(pgont_quantum, pgonr_quantum);
    poly_union_quantum = union(pgont_quantum, pgonr_quantum);
    
    poly_intersect_classical = intersect(pgont_classical, pgonr_classical);
    poly_union_classical = union(pgont_classical, pgonr_classical);
    
    figure(fig_classical_2)
    hold all;
    plot(pgont_classical);
    plot(pgonr_classical);
    plot(poly_intersect_classical);
    fn_fig_classical = sprintf('%s/lindbladian_evals_classical_poly_N(%d)_p(%0.4f)', path_figures, N, p);
    oqs_save_fig(fig_classical_2, fn_fig_classical);
    
    figure(fig_quantum_2)
    hold all;
    plot(pgont_quantum);
    plot(pgonr_quantum);
    plot(poly_intersect_quantum);
    fn_fig_quantum = sprintf('%s/lindbladian_evals_quantum_poly_N(%d)_p(%0.4f)', path_figures, N, p);
    oqs_save_fig(fig_quantum_2, fn_fig_quantum);
    
    iou_quantum(p_id) = area(poly_intersect_quantum) / area(poly_union_quantum);
    iou_classical(p_id) = area(poly_intersect_classical) / area(poly_union_classical);
    
end

close all