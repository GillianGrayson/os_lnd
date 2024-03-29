clear all;

N = 100;
num_seeds = 20000;

imag_lim = 1.3e-1;

pdf2d.x_num_bins = 201;
pdf2d.y_num_bins = 201;
pdf2d.x_label = '$Re(\chi'')$';
pdf2d.y_label = '$Im(\chi'')$';
pdf2d.x_bin_s = -300;
pdf2d.x_bin_f = 300;
pdf2d.y_bin_s = -300;
pdf2d.y_bin_f = 300;
pdf2d = oqs_pdf_2d_setup(pdf2d);

zs_all = zeros(N * num_seeds, 1);
s_id = 0;
f_id = 0;
for seed = 1:num_seeds
    tic
    
    rng(seed)
    
    fprintf('seed = %d\n', seed);
    
    G = Ginibre(seed, N);
    
    evals = eig(G);
    
    data2d = horzcat(real(evals), imag(evals));
    pdf2d = oqs_pdf_2d_update(pdf2d, data2d);
     
    evals_filtered = evals(abs(imag(evals)) >= imag_lim);
    num_passed_evals = size(evals_filtered, 1);
    fprintf('num_passed_evals = %d\n', num_passed_evals);
    
    s_id = f_id + 1;
    f_id = s_id + num_passed_evals - 1;
    
    zs = zeros(num_passed_evals, 1);
    neighbours_2D = horzcat(real(evals_filtered), imag(evals_filtered));
    for z_id = 1 : num_passed_evals
        target_2D = horzcat(real(evals_filtered(z_id)), imag(evals_filtered(z_id)));
        target = evals_filtered(z_id);
        distances = sqrt((neighbours_2D(:, 1) - target_2D(:, 1)).^2 + (neighbours_2D(:, 2) - target_2D(:, 2)).^2);
        [min_distances, order] = mink(distances, 3);
        nn = evals_filtered(order(2));
        nnn = evals_filtered(order(3));
        zs(z_id) = (nn - target) / (nnn - target);
    end
    
    num_zero_zs = sum(abs(zs) < 1e-4);
    fprintf('num_zero_zs = %d\n', num_zero_zs);
    
    zs_all(s_id : f_id) = zs;
    
    toc
end
zs_all = zs_all(1 : f_id);

total_num_passed_evals = f_id

pdf2d = oqs_pdf_2d_release(pdf2d);
fig = oqs_pdf_2d_plot(pdf2d);
fn_fig = sprintf('evals_N(%d)_num_seeds(%d)', N, num_seeds);
oqs_save_fig(fig, fn_fig);

pdf2dzs.x_num_bins = 301;
pdf2dzs.y_num_bins = 301;
pdf2dzs.x_label = '$Re(z)$';
pdf2dzs.y_label = '$Im(z)$';

pdf2dzs.x_bin_s = -1;
pdf2dzs.x_bin_f = 1;
pdf2dzs.y_bin_s = -1;
pdf2dzs.y_bin_f = 1;
pdf2dzs = oqs_pdf_2d_setup(pdf2dzs);
data2d = horzcat(real(zs_all), imag(zs_all));
pdf2dzs = oqs_pdf_2d_update(pdf2dzs, data2d);
pdf2dzs = oqs_pdf_2d_release(pdf2dzs);
fig = oqs_pdf_2d_plot(pdf2dzs);
fn_fig = sprintf('zs_N(%d)_num_seeds(%d)', N, num_seeds);
oqs_save_fig(fig, fn_fig)
