clear all;
path = '../../../source/cpp/os_lnd/os_lnd';
fig_path = 'E:/YandexDisk/Work/os_lnd/figures/xxz/local';

fontSize = 30;

num_spins = 7;

seed = 1;

Delta = 1.0;
W = 0.5;

mu = 0.001;

drv_type = 0;

ampl = 0.0;
freq = 6.2831853071795864769;
phase = 0.0;

j_id = 0;

num_states = 2.^num_spins;

suffix = sprintf('ns(%d)_prm(%0.4f_%0.4f_%0.4f_%d_%0.4f_%0.4f_%0.4f)', ...
    num_spins, ...
    Delta, ...
    W, ...
    mu, ...
    drv_type, ...
    ampl, ...
    freq, ...
    phase);

fn = sprintf('%s/rho_mtx_%s_seed(%d).txt', path, suffix, seed);
mtx_data = importdata(fn);

rho = zeros(num_states);
for st_id_1 = 1:num_states
    for st_id_2 = 1:num_states
        str_id = (st_id_1 - 1) * num_states + st_id_2;
        str = string(mtx_data(str_id));
        curr_data = sscanf(str, '(%e,%e)', 2);
        rho(st_id_1, st_id_2) = curr_data(1) + 1i * curr_data(2);
    end
end

xs = linspace(1, num_states, num_states);
abs_rho = abs(rho);

fig = figure;
propertyeditor('on');
imagesc(xs, xs, log10(abs_rho' + 1e-16));
set(gca, 'FontSize', fontSize);
xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(h, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', fontSize);
title(h, '$\log_{10}|\rho|$', 'FontSize', fontSize, 'interpreter','latex');
set(gca, 'TickLabelInterpreter', 'latex')
set(gca,'YDir','normal');
hold all;
fn = sprintf('%s/rho_mtx_%s_seed(%d)', fig_path, suffix, seed);
oqs_save_fig(gcf, fn)

fig_diag = figure;
propertyeditor('on');

colors = distinguishable_colors(num_spins);

for s_id = 1 : num_spins - 1
    fn = sprintf('%s/znd_mtx_%d_%s_seed(%d).txt', path, s_id - 1, suffix, seed);
    znd = oqs_read_sparse_mtx(fn, num_states);
    
    abs_znd = abs(znd);
 
    fig = figure;
    propertyeditor('on');
    imagesc(xs, xs, log10(abs_znd' + 1e-16));
    set(gca, 'FontSize', fontSize);
    xlabel('x', 'Interpreter', 'latex');
    ylabel('y', 'Interpreter', 'latex');
    colormap hot;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    set(gca, 'FontSize', fontSize);
    title(h, sprintf('$\\log_{10}|j_%d|$', s_id - 1), 'FontSize', fontSize, 'interpreter','latex');
    set(gca, 'TickLabelInterpreter', 'latex')
    set(gca,'YDir','normal');
    hold all;
    fn = sprintf('%s/znd_mtx_%d_%s_seed(%d)', fig_path, s_id - 1, suffix, seed);
    oqs_save_fig(gcf, fn);
    
    res_mtx = znd * rho;
    res_tr = trace(res_mtx);
    fprintf('j_id = %d trace = %0.16e + %0.16e i\n', s_id - 1, real(res_tr), imag(res_tr));
    fig = figure;
    propertyeditor('on');
    imagesc(xs, xs, log10(abs(res_mtx)' + 1e-16));
    set(gca, 'FontSize', fontSize);
    xlabel('x', 'Interpreter', 'latex');
    ylabel('y', 'Interpreter', 'latex');
    colormap hot;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    set(gca, 'FontSize', fontSize);
    title(h, sprintf('$\\log_{10}|R_%d|$', s_id - 1), 'FontSize', fontSize, 'interpreter','latex');
    set(gca, 'TickLabelInterpreter', 'latex')
    set(gca,'YDir','normal');
    hold all;
    fn = sprintf('%s/res_mtx_%d_%s_seed(%d)', fig_path, s_id - 1, suffix, seed);
    oqs_save_fig(gcf, fn);
    
    figure(fig_diag);
    h = scatter(xs, log10(real(diag(res_mtx))), 100, 'o', 'LineWidth',  1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', colors(s_id, :), 'MarkerEdgeAlpha', 0.6, 'MarkerFaceAlpha', 0.6);
    %h = plot(xs, log10(real(diag(res_mtx))), 'LineWidth', 1, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'white', 'MarkerSize', 8);
    legend(h, sprintf('k = %d, j = %0.2e', s_id - 1, real(res_tr)));
    set(gca, 'FontSize', fontSize);
    xlabel('x', 'Interpreter', 'latex');
    ylabel('$\log_{10}R_{n,n}$', 'Interpreter', 'latex');
    hold all;
    
end

figure(fig_diag);
legend(gca,'off');
legend('Location', 'SouthWest', 'NumColumns', 1, 'Interpreter', 'latex');
legend('FontSize', 20);
set(gca, 'TickLabelInterpreter', 'latex')
box on;
xlim([min(xs) - 0.5, max(xs) + 0.5])
fn = sprintf('%s/res_diag_%s_seed(%d)', fig_path, suffix, seed);
oqs_save_fig(gcf, fn);
