clear all;
addpath('../../../source/matlab/lib')

path = '/data/condmat/ivanchen/yusipov/os_lnd/mbl';
figures_path = '/home/ivanchen/yusipov/os_lnd/figures/mbl';

ns = 8;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

%Ws = linspace(1, 20, 20)';
Ws = [1, 20]';
U = 1.0;
J = 1.0;

seeds = linspace(1, 100, 100)';

num_Ws = size(Ws, 1);
num_seeds = size(seeds, 1);

x_num_bins = 101;
y_num_bins = 101;
x_label = '$Re(z)$';
y_label = '$Im(z)$';

N2 = (nchoosek(ns, ns/2)).^2

for W_id = 1:num_Ws
    
    W = Ws(W_id);
    fprintf('W = %0.16e\n', W);
    
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
    
    zs_all = zeros((N2 - 1) * size(seeds, 1), 1);
    rs_all = zeros((N2 - 1) * size(seeds, 1), 1);
    abses_all = zeros((N2 - 1) * size(seeds, 1), 1);
    angles_all = zeros((N2 - 1) * size(seeds, 1), 1);
    
    for seed_id = 1:num_seeds
	
        seed = seeds(seed_id);
		fprintf('seed = %d\n', seed);
		
        suffix = sprintf('ns(%d)_seed(%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)', ...
			ns, ...
            seed, ...
            diss_type, ...
            diss_phase, ...
            diss_gamma, ...
            W, ...
            U, ...
            J);
        
        fn = sprintf('%s/eigen_dense/ns_%d/diss_%d_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/seed_%d/lindbladian_evals_%s.txt', ...
            path, ...
            ns, ...
            diss_type, ...
            diss_phase, ...
            diss_gamma, ...
            W, ...
            U, ...
            J, ...
            seed, ...
            suffix);
        
        evals = zeros(N2, 1);
        data = importdata(fn);
        for str_id = 1:N2
            str = string(data(str_id));
            data2d = sscanf(str, '(%e,%e)', 2);
            evals(str_id) = data2d(1) + 1i * data2d(2);
        end
        [evals, order] = sort(evals, 'ComparisonMethod', 'abs');
        %fprintf('evals_first = %0.16e + i*%0.16e\n', real(evals(1)), imag(evals(1)));
        evals = evals(2:end);
        
        zs = zeros(N2 - 1, 1);
        rs = zeros(N2 - 1, 1);
        abses = zeros(N2 - 1, 1);
        angles = zeros(N2 - 1, 1);
        neighbours_2D = horzcat(real(evals), imag(evals));
        for z_id = 1 : N2 - 1
            target_2D = horzcat(real(evals(z_id)), imag(evals(z_id)));
            target = evals(z_id);
            distances = sqrt(sum(bsxfun(@minus, neighbours_2D, target_2D).^2, 2));
            [sorted_distances, order] = sort(distances);
            %fprintf('sorted_distances(1) = %0.16e\n', sorted_distances(1));
            nn = evals(order(2));
            nnn = evals(order(3));
            zs(z_id) = (nn - target) / (nnn - target);
            rs(z_id) = abs(zs(z_id)) * sign(angle(zs(z_id)));
            abses(z_id) = abs(zs(z_id));
            angles(z_id) = angle(zs(z_id));
        end
        
        
        s_id = (seed_id - 1) * (N2 - 1) + 1;
        f_id = seed_id * (N2 - 1);
        zs_all(s_id : f_id) = zs;
        rs_all(s_id : f_id) = rs;
        abses_all(s_id : f_id) = abses;
        angles_all(s_id : f_id) = angles;
    end
    
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

