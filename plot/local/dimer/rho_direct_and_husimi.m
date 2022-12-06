clear all;
path = '../../../source/cpp/os_lnd/os_lnd';
fig_path = 'E:/YandexDisk/Work/os_lnd/figures/dimer/local';

fontSize = 30;

task = 'odeint';

num_particles = 50;

diss_type = 1;
diss_gamma = 0.1;

E = 0.0;
U = 0.125;
J = 1.0;

drv_type = 1;
drv_ampl = 3.4;
drv_freq = 1.0;
drv_phase = 0.0;

num_states = num_particles + 1;

suffix = sprintf('np(%d)_diss(%d_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)', ...
    num_particles, ...
    diss_type, ...
    diss_gamma, ...
    E, ...
    U, ...
    J, ...
    drv_type, ...
    drv_ampl, ...
    drv_freq, ...
    drv_phase);

fn = sprintf('%s/rho_mtx_%s.txt', path, suffix);
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

tr = trace(rho)

xs = linspace(1, num_states, num_states);

fig = figure;
propertyeditor('on');
imagesc(xs, xs, log10(abs(rho)' + 1e-16));
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
fn = sprintf('%s/rho_mtx_log_%s', fig_path, suffix);
oqs_save_fig(gcf, fn)

fig = figure;
propertyeditor('on');
imagesc(xs, xs, abs(rho));
set(gca, 'FontSize', fontSize);
xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(h, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', fontSize);
title(h, '$|\rho|$', 'FontSize', fontSize, 'interpreter','latex');
set(gca, 'TickLabelInterpreter', 'latex')
set(gca,'YDir','normal');
hold all;
fn = sprintf('%s/rho_mtx_%s', fig_path, suffix);
oqs_save_fig(gcf, fn)

phi_size = 100;
phi_begin = 0;
phi_end = 2*pi;
phis = linspace(phi_begin, phi_end, phi_size)';

nu_size = 100;
nu_begin = 0;
nu_end = pi;
nus = linspace(nu_begin, nu_end, nu_size)';

tic
hus = husimi(nus, phis, transpose(rho));
toc

fig = figure;
propertyeditor(fig);
hLine = imagesc(nus, phis, real(hus'));
set(gca, 'FontSize', fontSize);
xlabel('$\vartheta$', 'Interpreter', 'latex');
set(gca, 'FontSize', fontSize);
ylabel('$\varphi$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', fontSize);
title(h, 'H', 'FontSize', fontSize, 'Interpreter', 'latex');
set(gca,'YDir','normal');
hold all;
fn = sprintf('%s/rho_mtx_husimi_%s', fig_path, suffix);
oqs_save_fig(gcf, fn)

