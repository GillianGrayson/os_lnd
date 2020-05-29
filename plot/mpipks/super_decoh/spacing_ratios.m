clear all;
addpath('../../../source/matlab/lib')

ps = [1.0]';

name = 'G_evals';

method = 'simple';
reshufle_type = 1;
G_type = 0;
N = 100;
seeds = linspace(1, 100, 100)';

x_num_bins = 301;
y_num_bins = 301;
x_label = '$Re(z)$';
y_label = '$Im(z)$';

N2 = N * N;

for p_id = 1:size(ps, 1)
    
    p = ps(p_id);
    
    fprintf('p = %0.16e\n', p);
    
    pdf2d.x_num_bins = x_num_bins;
    pdf2d.y_num_bins = y_num_bins;
    pdf2d.x_label = x_label;
    pdf2d.y_label = y_label;
    
    pdfrs.x_num_bins = x_num_bins;
    pdfrs.x_label = '$r$';
    
    pdfabses.x_num_bins = x_num_bins;
    pdfabses.x_label = '$|z|$';
    
    pdfangles.x_num_bins = x_num_bins;
    pdfangles.x_label = '$arg(z)$';
    
	evals_all = zeros((N2 - 1) * size(seeds, 1), 1);
    zs_all = zeros((N2 - 1) * size(seeds, 1), 1);
    rs_all = zeros((N2 - 1) * size(seeds, 1), 1);
    abses_all = zeros((N2 - 1) * size(seeds, 1), 1);
    angles_all = zeros((N2 - 1) * size(seeds, 1), 1);
    
    path = '/data/condmat/ivanchen/yusipov/os_lnd/super_decoh/eigen_dense';
    figures_path = '/home/ivanchen/yusipov/os_lnd/figures/super_decoh';
    
    N2 = N * N;
    
    for seed_id = 1:size(seeds, 1)
        
        seed = seeds(seed_id);
        fprintf('seed = %d\n', seed);
        
        suffix = sprintf('reshuffle(%d)_G(%d)_N(%d)_p(%0.10f)_seed(%d)', reshufle_type, G_type, N, p, seed);

        evals = zeros(N2, 1);
        fn_cpp = sprintf('%s/method_%s/G_type_%d/reshuffle_type_%d/N_%d/p_%0.10f/seed_%d/%s_%s.txt', path, method, G_type, reshufle_type, N, p, seed, name, suffix);
        if ~isfile(fn_cpp)
            fprintf('Warning! seed = %d\n', seed);
        end
        cpp_data = importdata(fn_cpp);
        for str_id = 1:N2
            str = string(cpp_data(str_id));
            data2d = sscanf(str, '(%e,%e)', 2);
            evals(str_id) = data2d(1) + 1i * data2d(2);
        end
        
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
            %distances = sqrt(sum(bsxfun(@minus, neighbours_2D, target_2D).^2, 2));
			%distances = zeros(N2 - 1, 1);
			%for tmp_id = 1 : N2 - 1
			%	tmp_2D = horzcat(real(evals(tmp_id)), imag(evals(tmp_id)));
			%	distances(tmp_id) = norm(target_2D - tmp_2D);
			%end
			distances = sqrt((neighbours_2D(:, 1) - target_2D(:, 1)).^2 + (neighbours_2D(:, 2) - target_2D(:, 2)).^2);
            [sorted_distances, order] = sort(distances);
			%fprintf('target(1) = %0.16e + %0.16e * i\n', real(evals(z_id)), imag(evals(z_id)));
            %fprintf('sorted_distances(1) = %0.16e\n', sorted_distances(1));
			%fprintf('sorted_distances(2) = %0.16e\n', sorted_distances(2));
			%fprintf('sorted_distances(3) = %0.16e\n', sorted_distances(3));
            nn = evals(order(2));
            nnn = evals(order(3));
            zs(z_id) = (nn - target) / (nnn - target);
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
    
    suffix = sprintf('%s_reshuffle(%d)_N(%d)_p(%0.10f)_numSeeds(%d)', name, reshufle_type, N, p, size(seeds, 1));
	
	fn_txt = sprintf('%s/zs_%s.txt', figures_path, suffix);
	%save(fn_txt, 'zs_all', '-double', '-tab');
    %dlmwrite(sprintf('%s/contour_%s.txt', figures_path, suffix), contour);
	%writematrix(zs_all, fn_txt, 'Delimiter','tab')
	fn_dat = sprintf('%s/evals_%s.dat', figures_path, suffix);
	dlmwrite(fn_dat, horzcat(real(evals_all(:)), imag(evals_all(:))))
		
	fid = fopen(fn_txt,'wt');
	for ii = 1:size(zs_all,1)
		fprintf(fid,'%0.16e\t%0.16e\n', real(zs_all(ii)), imag(zs_all(ii)));
	end
	fclose(fid)
    
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