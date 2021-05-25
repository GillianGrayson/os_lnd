clear all;
addpath('../../../source/matlab/lib')

path = '/data/biophys/denysov/yusipov/os_lnd/regular/super_decoh/eigen_dense';
figures_path = '/home/denysov/yusipov/os_lnd/figures/super_decoh';

is_calc_bow = 0;
is_save_files = 0;

ps = [0.0]';

region_eps = 1e-3;

name = 'lindbladian_evals';

method = 'origin';
reshuffle_type = 0;
G_type = 0;
N = 100;
num_seeds = 500;

x_num_bins = 301;
y_num_bins = 301;
x_label = '$Re(z)$';
y_label = '$Im(z)$';

N2 = N * N;

phase = 0 : pi/50 : 2*pi;
radius = 1.0;
l_circle_x = radius * cos(phase) + (1 - region_eps);
r_circle_x = radius * cos(phase) + (1 + region_eps);
l_circle_y = radius * sin(phase);
r_circle_y = radius * sin(phase);

bowstring_x = [0.5 - region_eps, 0.5 - region_eps, 0.5 + region_eps, 0.5 + region_eps, 0.5 - region_eps]';
bowstring_y = [-1, 1, 1, -1, -1]';

arrow_x = [-1, -1, 0, 0, -1];
arrow_y = [-region_eps, region_eps, region_eps, -region_eps, -region_eps];

for p_id = 1:size(ps, 1)
    
    p = ps(p_id);
    
    fprintf('p = %0.16e\n', p);
    
    pdf2d.x_num_bins = x_num_bins;
    pdf2d.y_num_bins = y_num_bins;
    pdf2d.x_label = x_label;
    pdf2d.y_label = y_label;
    
    pdf2d_zs_bow.x_num_bins = x_num_bins;
    pdf2d_zs_bow.y_num_bins = y_num_bins;
    pdf2d_zs_bow.x_label = x_label;
    pdf2d_zs_bow.y_label = y_label;
	
	pdf2d_evals_bow.x_num_bins = x_num_bins;
    pdf2d_evals_bow.y_num_bins = y_num_bins;
    pdf2d_evals_bow.x_label = '$Re(\lambda)$';
    pdf2d_evals_bow.y_label = '$Im(\lambda)$';
    
    pdf2d_zs_bowstring.x_num_bins = x_num_bins;
    pdf2d_zs_bowstring.y_num_bins = y_num_bins;
    pdf2d_zs_bowstring.x_label = x_label;
    pdf2d_zs_bowstring.y_label = y_label;
	
	pdf2d_evals_bowstring.x_num_bins = x_num_bins;
    pdf2d_evals_bowstring.y_num_bins = y_num_bins;
    pdf2d_evals_bowstring.x_label = '$Re(\lambda)$';
    pdf2d_evals_bowstring.y_label = '$Im(\lambda)$';
    
    pdf2d_zs_arrow.x_num_bins = x_num_bins;
    pdf2d_zs_arrow.y_num_bins = y_num_bins;
    pdf2d_zs_arrow.x_label = x_label;
    pdf2d_zs_arrow.y_label = y_label;
	
	pdf2d_evals_arrow.x_num_bins = x_num_bins;
    pdf2d_evals_arrow.y_num_bins = y_num_bins;
    pdf2d_evals_arrow.x_label = '$Re(\lambda)$';
    pdf2d_evals_arrow.y_label = '$Im(\lambda)$';
    
    pdfrs.x_num_bins = x_num_bins;
    pdfrs.x_label = '$r$';
    
    pdfabses.x_num_bins = x_num_bins;
    pdfabses.x_label = '$|z|$';
    
    pdfangles.x_num_bins = x_num_bins;
    pdfangles.x_label = '$arg(z)$';
    
	evals_all = zeros((N2 - 1) * num_seeds, 1);
    zs_all = zeros((N2 - 1) * num_seeds, 1);
    rs_all = zeros((N2 - 1) * num_seeds, 1);
    abses_all = zeros((N2 - 1) * num_seeds, 1);
    angles_all = zeros((N2 - 1) * num_seeds, 1);
    
    zs_bow = [];
    evals_bow = [];
    
    zs_bowstring = [];
    evals_bowstring = [];
    
    zs_arrow = [];
    evals_arrow = [];
    
    
    prefix = sprintf('method_%s/Gt_%d_ad_0/rt_%d/N_%d/p_%0.10f', ...
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
        num_seeds);
            
    fn = sprintf('%s/lindbladian_evals_%s.txt', aggr_path, suffix);

    lind_evals_all = importdata(fn);
	
	num_non_real = sum(abs(lind_evals_all(:,2)) > 1e-10)
    
    curr_fin_id = 0;
    curr_fin_slice_id = 0;
    
    for seed_id = 1:num_seeds
        
        fprintf('seed = %d\n', seed_id);
        
        evals = lind_evals_all((seed_id - 1) * N2 + 1 : seed_id * N2, 1) + 1i * lind_evals_all((seed_id - 1) * N2 + 1 : seed_id * N2, 2);
        [evals, order] = sort(evals, 'ComparisonMethod', 'abs');
        evals = evals(2:end);
        
        zs = zeros(N2 - 1, 1);
        rs = zeros(N2 - 1, 1);
        abses = zeros(N2 - 1, 1);
        angles = zeros(N2 - 1, 1);
        neighbours_2D = horzcat(real(evals), imag(evals));
        for z_id = 1 : N2 - 1
            target_2D = horzcat(real(evals(z_id)), imag(evals(z_id)));
            target = evals(z_id);
			distances = sqrt((neighbours_2D(:, 1) - target_2D(:, 1)).^2 + (neighbours_2D(:, 2) - target_2D(:, 2)).^2);
            % [min_distances, order] = mink(distances, 3);
            [sorted_distances, order] = sort(distances);
			nn = evals(order(2));
            nnn = evals(order(3));
            zs(z_id) = (nn - target) / (nnn - target);
            
            if is_calc_bow == 1
                l_in = inpolygon(real(zs(z_id)), imag(zs(z_id)), l_circle_x, l_circle_y);
                r_in = inpolygon(real(zs(z_id)), imag(zs(z_id)), r_circle_x, r_circle_y);
                if (l_in == 1) && (r_in == 0)
                    zs_bow = vertcat(zs_bow, zs(z_id));
                    evals_bow = vertcat(evals_bow, [target; nn; nnn]);
                end

                bow_string_in = inpolygon(real(zs(z_id)), imag(zs(z_id)), bowstring_x, bowstring_y);
                if (bow_string_in == 1)
                    zs_bowstring = vertcat(zs_bowstring, zs(z_id));
                    evals_bowstring = vertcat(evals_bowstring, [target; nn; nnn]);
                end

                arrow_in = inpolygon(real(zs(z_id)), imag(zs(z_id)), arrow_x, arrow_y);
                if (arrow_in == 1)
                    zs_arrow = vertcat(zs_arrow, zs(z_id));
                    evals_arrow = vertcat(evals_arrow, [target; nn; nnn]);
                end
                
            end
            
            rs(z_id) = abs(zs(z_id)) * sign(angle(zs(z_id)));
            abses(z_id) = abs(zs(z_id));
            angles(z_id) = angle(zs(z_id));
        end
        
        s_id = (seed_id - 1) * (N2 - 1) + 1;
        f_id = seed_id * (N2 - 1);
		evals_all(s_id : f_id) = evals;
        zs_all(s_id : f_id) = zs;
        rs_all(s_id : f_id) = rs;
        abses_all(s_id : f_id) = abses;
        angles_all(s_id : f_id) = angles;
    end
    
    suffix = sprintf('%s_reshuffle(%d)_N(%d)_p(%0.10f)_numSeeds(%d)', name, reshuffle_type, N, p, num_seeds);
	
    if is_save_files == 1
        fn_txt = sprintf('%s/zs_%s.txt', figures_path, suffix);
        save(fn_txt, 'zs_all', '-double', '-tab');
        dlmwrite(sprintf('%s/contour_%s.txt', figures_path, suffix), contour);
        writematrix(zs_all, fn_txt, 'Delimiter','tab')
        fn_dat = sprintf('%s/evals_%s.dat', figures_path, suffix);
        dlmwrite(fn_dat, horzcat(real(evals_all(:)), imag(evals_all(:))))

        fid = fopen(fn_txt,'wt');
        for ii = 1:size(zs_all,1)
        	fprintf(fid,'%0.16e\t%0.16e\n', real(zs_all(ii)), imag(zs_all(ii)));
        end
        fclose(fid)
    end
    
    pdf2d.x_bin_s = min(real(zs_all));
    pdf2d.x_bin_f = max(real(zs_all));
    pdf2d.y_bin_s = min(imag(zs_all));
    pdf2d.y_bin_f = max(imag(zs_all));
    pdf2d = oqs_pdf_2d_setup(pdf2d);
    data2d = horzcat(real(zs_all), imag(zs_all));
    pdf2d = oqs_pdf_2d_update(pdf2d, data2d);
    pdf2d = oqs_pdf_2d_release(pdf2d);
    fig = oqs_pdf_2d_plot(pdf2d);
    fn_fig = sprintf('%s/zs_%s', figures_path, suffix);
    oqs_save_fig(fig, fn_fig)
	
	fig = figure;
	imagesc(pdf2d.x_bin_centers, pdf2d.y_bin_centers, log10(pdf2d.pdf' + 1e-16));
	set(gca, 'FontSize', 30);
	xlabel(pdf2d.x_label, 'Interpreter', 'latex');
	set(gca, 'FontSize', 30);
	ylabel(pdf2d.y_label, 'Interpreter', 'latex');
	colormap hot;
	h = colorbar;
	set(gca, 'FontSize', 30);
	title(h, '$\log_{10}PDF$', 'FontSize', 33, 'interpreter','latex');
	set(gca,'YDir','normal');
	hold all;

	fn_fig = sprintf('%s/zs_log_%s', figures_path, suffix);
	oqs_save_fig(fig, fn_fig)
    
    if is_save_files == 1
        fn_txt = sprintf('%s/zs_pdf_%s.txt', figures_path, suffix);
        fid = fopen(fn_txt,'wt');
        for x_id = 1:size(pdf2d.x_bin_centers, 2)
            for y_id = 1:size(pdf2d.y_bin_centers, 2)
                fprintf(fid,'%0.16e\t%0.16e\t%0.16e\n', pdf2d.x_bin_centers(x_id), pdf2d.y_bin_centers(y_id), pdf2d.pdf(x_id, y_id));
            end
        end
        fclose(fid)
    end
    
    if is_calc_bow == 1
        pdf2d_zs_bow.x_bin_s = min(real(zs_all));
        pdf2d_zs_bow.x_bin_f = max(real(zs_all));
        pdf2d_zs_bow.y_bin_s = min(imag(zs_all));
        pdf2d_zs_bow.y_bin_f = max(imag(zs_all));
        pdf2d_zs_bow = oqs_pdf_2d_setup(pdf2d_zs_bow);
        data2d = horzcat(real(zs_bow), imag(zs_bow));
        pdf2d_zs_bow = oqs_pdf_2d_update(pdf2d_zs_bow, data2d);
        pdf2d_zs_bow.pdf = pdf2d_zs_bow.pdf / (pdf2d.inc_count * pdf2d.x_bin_shift * pdf2d.y_bin_shift);
        pdf2d_zs_bow.norm = sum(sum(pdf2d_zs_bow.pdf)) * pdf2d_zs_bow.x_bin_shift * pdf2d_zs_bow.y_bin_shift;
        fprintf('pdf_norm = %0.16e\n', pdf2d_zs_bow.norm);
        fig = oqs_pdf_2d_plot(pdf2d_zs_bow);
        fn_fig = sprintf('%s/zs_bow_%s', figures_path, suffix);
        oqs_save_fig(fig, fn_fig)

        pdf2d_evals_bow.x_bin_s = min(real(evals_bow));
        pdf2d_evals_bow.x_bin_f = max(real(evals_bow));
        pdf2d_evals_bow.y_bin_s = min(imag(evals_bow));
        pdf2d_evals_bow.y_bin_f = max(imag(evals_bow));
        pdf2d_evals_bow = oqs_pdf_2d_setup(pdf2d_evals_bow);
        data2d = horzcat(real(evals_bow), imag(evals_bow));
        pdf2d_evals_bow = oqs_pdf_2d_update(pdf2d_evals_bow, data2d);
        pdf2d_evals_bow = oqs_pdf_2d_release(pdf2d_evals_bow);
        fig = oqs_pdf_2d_plot(pdf2d_evals_bow);
        fn_fig = sprintf('%s/evals_bow_%s', figures_path, suffix);
        oqs_save_fig(fig, fn_fig)

        pdf2d_zs_bowstring.x_bin_s = min(real(zs_all));
        pdf2d_zs_bowstring.x_bin_f = max(real(zs_all));
        pdf2d_zs_bowstring.y_bin_s = min(imag(zs_all));
        pdf2d_zs_bowstring.y_bin_f = max(imag(zs_all));
        pdf2d_zs_bowstring = oqs_pdf_2d_setup(pdf2d_zs_bowstring);
        data2d = horzcat(real(zs_bowstring), imag(zs_bowstring));
        pdf2d_zs_bowstring = oqs_pdf_2d_update(pdf2d_zs_bowstring, data2d);
        pdf2d_zs_bowstring.pdf = pdf2d_zs_bowstring.pdf / (pdf2d.inc_count * pdf2d.x_bin_shift * pdf2d.y_bin_shift);
        pdf2d_zs_bowstring.norm = sum(sum(pdf2d_zs_bowstring.pdf)) * pdf2d_zs_bowstring.x_bin_shift * pdf2d_zs_bowstring.y_bin_shift;
        fprintf('pdf_norm = %0.16e\n', pdf2d_zs_bowstring.norm);
        fig = oqs_pdf_2d_plot(pdf2d_zs_bowstring);
        fn_fig = sprintf('%s/zs_bowstring_%s', figures_path, suffix);
        oqs_save_fig(fig, fn_fig)

        pdf2d_evals_bowstring.x_bin_s = min(real(evals_bow));
        pdf2d_evals_bowstring.x_bin_f = max(real(evals_bow));
        pdf2d_evals_bowstring.y_bin_s = min(imag(evals_bow));
        pdf2d_evals_bowstring.y_bin_f = max(imag(evals_bow));
        pdf2d_evals_bowstring = oqs_pdf_2d_setup(pdf2d_evals_bowstring);
        data2d = horzcat(real(evals_bowstring), imag(evals_bowstring));
        pdf2d_evals_bowstring = oqs_pdf_2d_update(pdf2d_evals_bowstring, data2d);
        pdf2d_evals_bowstring = oqs_pdf_2d_release(pdf2d_evals_bowstring);
        fig = oqs_pdf_2d_plot(pdf2d_evals_bowstring);
        fn_fig = sprintf('%s/evals_bowstring_%s', figures_path, suffix);
        oqs_save_fig(fig, fn_fig)

        pdf2d_zs_arrow.x_bin_s = min(real(zs_all));
        pdf2d_zs_arrow.x_bin_f = max(real(zs_all));
        pdf2d_zs_arrow.y_bin_s = min(imag(zs_all));
        pdf2d_zs_arrow.y_bin_f = max(imag(zs_all));
        pdf2d_zs_arrow = oqs_pdf_2d_setup(pdf2d_zs_arrow);
        data2d = horzcat(real(zs_arrow), imag(zs_arrow));
        pdf2d_zs_arrow = oqs_pdf_2d_update(pdf2d_zs_arrow, data2d);
        pdf2d_zs_arrow.pdf = pdf2d_zs_arrow.pdf / (pdf2d.inc_count * pdf2d.x_bin_shift * pdf2d.y_bin_shift);
        pdf2d_zs_arrow.norm = sum(sum(pdf2d_zs_arrow.pdf)) * pdf2d_zs_arrow.x_bin_shift * pdf2d_zs_arrow.y_bin_shift;
        fprintf('pdf_norm = %0.16e\n', pdf2d_zs_arrow.norm);
        fig = oqs_pdf_2d_plot(pdf2d_zs_arrow);
        fn_fig = sprintf('%s/zs_arrow_%s', figures_path, suffix);
        oqs_save_fig(fig, fn_fig)

        pdf2d_evals_arrow.x_bin_s = min(real(evals_bow));
        pdf2d_evals_arrow.x_bin_f = max(real(evals_bow));
        pdf2d_evals_arrow.y_bin_s = min(imag(evals_bow));
        pdf2d_evals_arrow.y_bin_f = max(imag(evals_bow));
        pdf2d_evals_arrow = oqs_pdf_2d_setup(pdf2d_evals_arrow);
        data2d = horzcat(real(evals_arrow), imag(evals_arrow));
        pdf2d_evals_arrow = oqs_pdf_2d_update(pdf2d_evals_arrow, data2d);
        pdf2d_evals_arrow = oqs_pdf_2d_release(pdf2d_evals_arrow);
        fig = oqs_pdf_2d_plot(pdf2d_evals_arrow);
        fn_fig = sprintf('%s/evals_arrow_%s', figures_path, suffix);
        oqs_save_fig(fig, fn_fig)
        
    end
    
    pdfrs.x_bin_s = min(rs_all);
    pdfrs.x_bin_f = max(rs_all);
    pdfrs = oqs_pdf_1d_setup(pdfrs);
    pdfrs = oqs_pdf_1d_update(pdfrs, rs_all);
    pdfrs = oqs_pdf_1d_release(pdfrs);
    fig = oqs_pdf_1d_plot(pdfrs);
    fn_fig = sprintf('%s/rs_%s', figures_path, suffix);
    oqs_save_fig(fig, fn_fig);
    
    pdfabses.x_bin_s = min(abses_all);
    pdfabses.x_bin_f = max(abses_all);
    pdfabses = oqs_pdf_1d_setup(pdfabses);
    pdfabses = oqs_pdf_1d_update(pdfabses, abses_all);
    pdfabses = oqs_pdf_1d_release(pdfabses);
    fig = oqs_pdf_1d_plot(pdfabses);
    fn_fig = sprintf('%s/abses_%s', figures_path, suffix);
    oqs_save_fig(fig, fn_fig);
    
    pdfangles.x_bin_s = min(angles_all);
    pdfangles.x_bin_f = max(angles_all);
    pdfangles = oqs_pdf_1d_setup(pdfangles);
    pdfangles = oqs_pdf_1d_update(pdfangles, angles_all);
    pdfangles = oqs_pdf_1d_release(pdfangles);
    fig = oqs_pdf_1d_plot(pdfangles);
    fn_fig = sprintf('%s/angles_%s', figures_path, suffix);
    oqs_save_fig(fig, fn_fig);
    
end