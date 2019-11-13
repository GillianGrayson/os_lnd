clear all;

data_path = '/data/biophys/denysov/yusipov/os_lnd/dimer';
figures_path = '/home/denysov/yusipov/os_lnd/figures/dimer';

task = 'lindbladian_odeint_rk4';

num_particles = 100;

diss_type = 1;
diss_gamma = 0.1;

E = 0.0;
Us = [0.01 : 0.01 : 1.0]';
J = 1.0;

drv_type = 1;
drv_ampl = 3.4;
drv_freq = 1.0;
drv_phase = 0.0;

num_states = num_particles + 1;
states = [0:1:num_particles]';
diag_rho = zeros(size(Us, 1), num_particles + 1);

for U_id = 1 : size(Us, 1)
    
    U = Us(U_id)
    
    local_path = sprintf('%s/np_%d/diss_%d_%0.4f/prm_%0.4f_%0.4f_%0.4f/drv_%d_%0.4f_%0.4f_%0.4f', ...
        task, ...
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
    
    fn = sprintf('%s/%s/rho_mtx_%s.txt', data_path, local_path, suffix);
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
    
    curr_diag_rho = abs(diag(rho));
    curr_diag_rho = curr_diag_rho / max(curr_diag_rho);
    
    diag_rho(U_id, :) = curr_diag_rho;
end

fig = figure;
hLine = imagesc(states, Us, diag_rho');
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$PDF$','Interpreter', 'latex');
set(gca,'YDir','normal');

suffix = sprintf('np(%d)_diss(%d_%0.4f)_prm(%0.4f_var_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)', ...
    num_particles, ...
    diss_type, ...
    diss_gamma, ...
    E, ...
    J, ...
    drv_type, ...
    drv_ampl, ...
    drv_freq, ...
    drv_phase);

savefig(sprintf('%s/%s/diag_rho_from_U_%s.fig', figures_path, task, suffix));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/%s/diag_rho_from_U_%s.pdf', figures_path, task, suffix));