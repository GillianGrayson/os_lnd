clear all;

N = 2000;
num_seeds = 10;

path_border_line = 'E:/YandexDisk/Work/os_lnd/data';
theory_data_classical = importdata(sprintf('%s/borderline_classical.dat', path_border_line));

pdf2d.x_num_bins = 201;
pdf2d.y_num_bins = 201;
pdf2d.x_label = 'Re$(\lambda)$';
pdf2d.y_label = 'Im$(\lambda)$';
pdf2d.x_bin_s = -4;
pdf2d.x_bin_f = 4;
pdf2d.y_bin_s = -1;
pdf2d.y_bin_f = 1;
pdf2d = oqs_pdf_2d_setup(pdf2d);
all_evals = zeros(N, num_seeds);
for seed = 1:num_seeds
    tic
    
    fprintf('seed = %d\n', seed);
    
    A = zeros(N);
    G = getG(N, seed);
    for i = 1:N
        for j = 1:N
            A(i,j) = (abs(G(i,j)))^2;
        end
    end
    K = A - eye(N).*sum(A,2);
    evals = eig(K);
    all_evals(:, seed) = evals;
    
    evals = evals(2:end);
    evals = sqrt(N) * (real(evals) + 1) + 1i * sqrt(N) * imag(evals);
    
    data2d = horzcat(real(evals), imag(evals));
    pdf2d = oqs_pdf_2d_update(pdf2d, data2d);

    toc
end

fig = oqs_pdf_2d_plot(pdf2d);
hold all;
xt1 = theory_data_classical(:, 1);
yt1 = theory_data_classical(:, 2);
xt2 = flip(xt1);
yt2 = -flip(yt1);
xt = vertcat(xt1, xt2);
yt = vertcat(yt1, yt2);
hLine = plot(xt, yt, 'LineWidth', 3, 'Color', 'cyan');
fn_fig = sprintf('evals');
oqs_save_fig(fig, fn_fig);

