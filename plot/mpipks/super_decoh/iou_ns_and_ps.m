clear all;
addpath('../../../source/matlab/lib')

task = 'eigen_dense';
path = sprintf('/data/condmat/ivanchen/yusipov/os_lnd/regular/super_decoh/%s', task);

figures_path = '/home/ivanchen/yusipov/os_lnd/figures/super_decoh';

num_hits = 0;

method = 'simple';
G_type = 0;
reshuffle_type = 0;

ps = logspace(-3, 0, 31)';
%Ns = floor(logspace(1, 2, 11)');
Ns = vertcat(floor(logspace(1, 2, 11)'), 200);
total_num_evals = 100000;
num_seeds = zeros(size(Ns, 1), 1);

point_x = zeros(size(Ns, 1), 1);
point_y = zeros(size(Ns, 1), 1);

fig_global = figure;

colors = distinguishable_colors(size(Ns, 1));

for N_id = 1:size(Ns, 1)

	color = colors(N_id, :)
    
    N = Ns(N_id);
    N2 = N * N;
    
	if N == 200
		num_seeds(N_id) = 10;
	else
		num_seeds(N_id) = ceil(total_num_evals / (Ns(N_id).^2));
    end
	
    fprintf('N = %d\n', N);
    fprintf('num_seeds = %d\n', num_seeds(N_id));
    
    theory_data_classical = importdata(sprintf('borderline_classical.dat'));
    theory_data_quantum = importdata(sprintf('borderline_quantum.dat'));
    
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

        all_evals_classical = zeros((N2 - 1) * num_seeds(N_id), 1);
        all_evals_quantum = zeros((N2 - 1) * num_seeds(N_id), 1);
  
        prefix = sprintf('method_%s/G_type_%d_ad_0/reshuffle_type_%d/N_%d/p_%0.10f', ...
            method, ...
            G_type, ...
            reshuffle_type, ...
            N, ...
            p);
        
        aggr_path = sprintf('%s/aggregator/%s', ...
            path, ...
            prefix);
        mkdir(aggr_path);
        
        suffix = sprintf('reshuffle(%d)_G(%d)_ad(0)_N(%d)_p(%0.10f)_seeds(%d_%d_%d)', ...
            reshuffle_type, ...
            G_type, ...
            N, ...
            p, ...
            1, ...
            1, ...
            num_seeds(N_id));
        
        fn_txt = sprintf('%s/lindbladian_evals_%s.txt', aggr_path, suffix);
        lind_evals_all = importdata(fn_txt);
        
        for seed_id = 1:num_seeds(N_id)
            
            seed = seed_id;
            
            evals = lind_evals_all((seed_id - 1) * N2 + 1 : seed_id * N2, 1)  + 1i * lind_evals_all((seed_id - 1) * N2 + 1 : seed_id * N2, 2);
            
            [evals, order] = sort(evals, 'ComparisonMethod', 'abs');
            evals = evals(2:end);
            evals_classical = N * sqrt(N) * (real(evals) + 1) + 1i * N * sqrt(N) * imag(evals);
            evals_quantum = N / p * (real(evals) + 1) + 1i * N / p * imag(evals);

            s_id = (seed_id - 1) * (N2 - 1) + 1;
            f_id = seed_id * (N2 - 1);
            all_evals_classical(s_id : f_id) = evals_classical;
            all_evals_quantum(s_id : f_id) = evals_quantum;
            
        end
        
        suffix_classical = sprintf('numHits(%d)_type(classical)_reshuffle(%d)_N(%d)_p(%0.10f)_numSeeds(%d)', num_hits, reshuffle_type, N, p, num_seeds(N_id));
        suffix_quantum = sprintf('numHits(%d)_type(quantum)_reshuffle(%d)_N(%d)_p(%0.10f)_numSeeds(%d)', num_hits, reshuffle_type, N, p, num_seeds(N_id));
        
        pdf2d_classical.x_bin_s = min(real(all_evals_classical)) - 1e-16;
        pdf2d_classical.x_bin_f = max(real(all_evals_classical)) + 1e-16;
        pdf2d_classical.y_bin_s = min(imag(all_evals_classical)) - 1e-16;
        pdf2d_classical.y_bin_f = max(imag(all_evals_classical)) + 1e-16;
        pdf2d_classical = oqs_pdf_2d_setup(pdf2d_classical);
        data2d = horzcat(real(all_evals_classical), imag(all_evals_classical));
        pdf2d_classical = oqs_pdf_2d_update(pdf2d_classical, data2d);
        fprintf('num_evals = %d\n', pdf2d_classical.inc_count);
        pdf2d_classical = oqs_pdf_2d_release(pdf2d_classical);
        fig_classical = oqs_pdf_2d_plot(pdf2d_classical);
        fn_fig_classical = sprintf('%s/lindbladian_evals_%s', figures_path, suffix_classical);
        
        pdf2d_quantum.x_bin_s = min(real(all_evals_quantum)) - 1e-16;
        pdf2d_quantum.x_bin_f = max(real(all_evals_quantum)) + 1e-16;
        pdf2d_quantum.y_bin_s = min(imag(all_evals_quantum)) - 1e-16;
        pdf2d_quantum.y_bin_f = max(imag(all_evals_quantum)) + 1e-16;
        pdf2d_quantum = oqs_pdf_2d_setup(pdf2d_quantum);
        data2d = horzcat(real(all_evals_quantum), imag(all_evals_quantum));
        pdf2d_quantum = oqs_pdf_2d_update(pdf2d_quantum, data2d);
        fprintf('num_evals = %d\n', pdf2d_quantum.inc_count);
        pdf2d_quantum = oqs_pdf_2d_release(pdf2d_quantum);
        fig_quantum = oqs_pdf_2d_plot(pdf2d_quantum);
        fn_fig_quantum = sprintf('%s/lindbladian_evals_%s', figures_path, suffix_quantum);
        
        xt1 = theory_data_quantum(:, 1);
        yt1 = theory_data_quantum(:, 2);
        xt2 = flip(xt1);
        yt2 = -flip(yt1);
        xt = vertcat(xt1, xt2);
        yt = vertcat(yt1, yt2);
        pgont_quantum = polyshape(xt, yt);
        
        xt1 = theory_data_classical(:, 1);
        yt1 = theory_data_classical(:, 2);
        xt2 = flip(xt1);
        yt2 = -flip(yt1);
        xt = vertcat(xt1, xt2);
        yt = vertcat(yt1, yt2);
        pgont_classical = polyshape(xt, yt);
        
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
        
        figure(fig_quantum)
        hold all;
        plot(pgont_quantum);
        plot(pgonr_quantum);
        plot(poly_intersect_quantum);
        
        figure(fig_classical)
        hold all;
        plot(pgont_classical);
        plot(pgonr_classical);
        plot(poly_intersect_classical);
        
        iou_quantum(p_id) = area(poly_intersect_quantum) / area(poly_union_quantum);
        iou_classical(p_id) = area(poly_intersect_classical) / area(poly_union_classical);
        
        fprintf('iou_area_quantum = %0.16e\n', iou_quantum(p_id));
        fprintf('iou_area_lassical = %0.16e\n', iou_classical(p_id));
        
        %oqs_save_fig(fig_quantum, fn_fig_quantum);
        %oqs_save_fig(fig_classical, fn_fig_classical);
    end
    
    figure(fig_global)
    hline = plot(ps, iou_quantum, 'LineWidth', 2, 'Color', color);
    color = get(hline,'Color');
    legend(hline, sprintf('N=%d', N))
    set(gca, 'FontSize', 30);
    xlabel('p', 'Interpreter', 'latex');
    set(gca, 'FontSize', 30);
    ylabel('IoU', 'Interpreter', 'latex');
    set(gca,'XScale','log');
    %set(gca,'YScale','log');
	hold all;

    %mean_iou = (max(iou_quantum) + min(iou_quantum)) * 0.5;
    %[xi, yi] = polyxpoly(ps, iou_quantum, [ps(1) ps(end)], [mean_iou mean_iou]);
    [xi, yi] = polyxpoly(ps, iou_quantum, ps, iou_classical);
	if size(xi, 1) > 1
		xi = xi(end);
		yi = yi(end);
	end
    hline = plot(xi, yi, 'o', 'MarkerSize', 10, 'MarkerFaceColor', color);
    hline.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold all;
	
	[val, idx]=min(abs(ps - xi));
	start_index = idx;
	if ps(idx) < xi
		start_index = idx + 1;
	end	
    %hline = plot(vertcat(xi, ps(start_index:end)), vertcat(yi, iou_classical(start_index:end)), 'LineWidth', 2, 'Color', color);
    %hline.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %hold all;
	
	hline = plot(ps(1:end), iou_classical(1:end), 'LineWidth', 2, 'Color', color, 'LineStyle', '--');
    hline.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold all;
    
    point_x(N_id) = xi;
    point_y(N_id) = yi;
end

suffix_quantum = sprintf('numHits(%d)_reshuffle(%d)_N(var)_p(var)', num_hits, reshuffle_type);
fn_fig_quantum = sprintf('%s/iou_%s', figures_path, suffix_quantum);
oqs_save_fig(fig_global, fn_fig_quantum)

fig_quantum = figure;
hline = plot(Ns, point_x, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('N', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$p_{trans}$', 'Interpreter', 'latex');
set(gca,'XScale','log');
set(gca,'YScale','log');
fn_fig_quantum = sprintf('%s/point_%s', figures_path, suffix_quantum);
oqs_save_fig(fig_quantum, fn_fig_quantum)
