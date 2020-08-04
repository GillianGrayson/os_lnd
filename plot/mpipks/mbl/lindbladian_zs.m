clear all;
addpath('../../../source/matlab/lib')

task = 'eigen_dense';

path = sprintf('/data/condmat/ivanchen/yusipov/os_lnd/mbl/%s', task);

figures_path = '/home/ivanchen/yusipov/os_lnd/figures/mbl';

without_line = 0;

ns = 8;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

Ws = [1, 20]';
U = 1.0;
J = 1.0;

global_seeds = 10000;
target_seeds = 1000;

seeds = linspace(1, target_seeds, target_seeds)';

num_Ws = size(Ws, 1);
num_seeds = size(seeds, 1);

x_num_bins = 201;
y_num_bins = 201;
x_label = '$Re(z)$';
y_label = '$Im(z)$';

N = (nchoosek(ns, ns/2));
N2 = N^2;

for W_id = 1:num_Ws
    
    W = Ws(W_id);
    fprintf('W = %0.16e\n', W);
    
    pdf2d.x_num_bins = x_num_bins;
    pdf2d.y_num_bins = y_num_bins;
    pdf2d.x_label = x_label;
    pdf2d.y_label = y_label;
    
    pdfrs.x_num_bins = x_num_bins;
    pdfrs.x_label = '$r$';
    
    pdfzs.x_num_bins = x_num_bins;
    pdfzs.x_label = '$z$';
    
    pdfabses.x_num_bins = x_num_bins;
    pdfabses.x_label = '$|z|$';
    
    pdfangles.x_num_bins = x_num_bins;
    pdfangles.x_label = '$arg(z)$';
    
    zs_all = zeros(N2 * size(seeds, 1), 1);
    rs_all = zeros(N2 * size(seeds, 1), 1);
    abses_all = zeros(N2 * size(seeds, 1), 1);
    angles_all = zeros(N2 * size(seeds, 1), 1);
    
    zs_slice_all = zeros(N2 * size(seeds, 1), 1);
       
    prefix = sprintf('ns_%d/diss_%d_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f', ...
        ns, ...
        diss_type, ...
        diss_phase, ...
        diss_gamma, ...
        W, ...
        U, ...
        J);
    
    aggr_path = sprintf('%s/aggregator/%s', ...
        path, ...
        prefix);
    
    suffix = sprintf('ns(%d)_seeds(%d_%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)', ...
        ns, ...
        1, ...
        1, ...
        global_seeds, ...
        diss_type, ...
        diss_phase, ...
        diss_gamma, ...
        W, ...
        U, ...
        J);
    
    fn = sprintf('%s/lindbladian_evals_%s.txt', aggr_path, suffix);

    lind_evals_all = importdata(fn);
    
    curr_fin_id = 0;
    curr_fin_slice_id = 0;
    
    for seed_id = 1:num_seeds
	
        seed = seeds(seed_id);
		fprintf('seed = %d\n', seed);
		
        evals = lind_evals_all((seed_id - 1) * N2 + 1 : seed_id * N2, 1) + 1i * lind_evals_all((seed_id - 1) * N2 + 1 : seed_id * N2, 2);
        [evals, order] = sort(evals, 'ComparisonMethod', 'abs');
        evals = evals(2:end);
        
        if without_line == 1
            evals = evals(abs(imag(evals(:))) > 1e-2); 
        end
        
        s_id = curr_fin_id + 1;
        f_id = curr_fin_id + size(evals, 1);
        curr_fin_id = f_id;
        
        zs = zeros(size(evals, 1), 1);
        rs = zeros(size(evals, 1), 1);
        abses = zeros(size(evals, 1), 1);
        angles = zeros(size(evals, 1), 1);
        
        zs_slice = [];
        
        neighbours_2D = horzcat(real(evals), imag(evals));
        for z_id = 1 : size(evals, 1)
            target_2D = horzcat(real(evals(z_id)), imag(evals(z_id)));
            target = evals(z_id);
            distances = sqrt(sum(bsxfun(@minus, neighbours_2D, target_2D).^2, 2));
            [sorted_distances, order] = sort(distances);
            nn = evals(order(2));
            nnn = evals(order(3));
            zs(z_id) = (nn - target) / (nnn - target);
            
            if abs(imag(zs(z_id))) < 1e-2
                zs_slice = vertcat(zs_slice, zs(z_id));
            end
            
            rs(z_id) = abs(zs(z_id)) * sign(angle(zs(z_id)));
            abses(z_id) = abs(zs(z_id));
            angles(z_id) = angle(zs(z_id));
        end
        
        zs_all(s_id : f_id) = zs;
        rs_all(s_id : f_id) = rs;
        abses_all(s_id : f_id) = abses;
        angles_all(s_id : f_id) = angles;
        
        zs_slice_all(curr_fin_slice_id + 1 : curr_fin_slice_id + size(zs_slice, 1)) = zs_slice;
        curr_fin_slice_id = curr_fin_slice_id + size(zs_slice, 1);
        
    end
    
	fprintf('curr_fin_id = %d\n', curr_fin_id);
	
    zs_all = zs_all(1 : curr_fin_id);
    rs_all = rs_all(1 : curr_fin_id);
    abses_all = abses_all(1 : curr_fin_id);
    angles_all = angles_all(1 : curr_fin_id);
    
    zs_slice_all = zs_slice_all(1 : curr_fin_slice_id);
    
    suffix = sprintf('ns(%d)_numSeeds(%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)', ...
        ns, ...
        num_seeds, ...
        diss_type, ...
        diss_phase, ...
        diss_gamma, ...
        W, ...
        U, ...
        J);
    
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
    
    pdfzs.x_bin_s = -1;
    pdfzs.x_bin_f = +1;
    pdfzs = oqs_pdf_1d_setup(pdfzs);
    pdfzs = oqs_pdf_1d_update(pdfzs, real(zs_slice_all));
    pdfzs = oqs_pdf_1d_release(pdfzs);
    fig = oqs_pdf_1d_plot(pdfzs);
    fn_fig = sprintf('%s/zs_slice_%s', figures_path, suffix);
    oqs_save_fig(fig, fn_fig);
    
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

