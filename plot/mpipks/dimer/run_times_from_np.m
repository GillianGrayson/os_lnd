clear all;

data_path = '/data/biophys/denysov/yusipov/os_lnd/dimer';
figures_path = '/home/denysov/yusipov/os_lnd/figures/dimer';

task = 'lindbladian_odeint_rk4';

start_observed_period = 0;
finish_observed_period = 1;
step = (2.0 * pi) * 0.0001;

nums_particles = [10:10:1000]';

diss_type = 1;
diss_gamma = 0.1;

E = 0.0;
U = 0.01;
J = 1.0;

drv_type = 1;
drv_ampl = 3.4;
drv_freq = 1.0;
drv_phase = 0.0;

init_times = zeros(size(nums_particles, 1), 1);
odeint_times = zeros(size(nums_particles, 1), 1);

for n_id = 1 : size(nums_particles, 1)
    
    num_particles = nums_particles(n_id)
    
    local_path = sprintf('%s_%d_%d_%0.2e/np_%d/diss_%d_%0.4f/prm_%0.4f_%0.4f_%0.4f/drv_%d_%0.4f_%0.4f_%0.4f', ...
        task, ...
		start_observed_period, ...
		finish_observed_period, ...
		step, ...
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
    
    fn = sprintf('%s/%s/run_times_%s.txt', data_path, local_path, suffix);
    data = importdata(fn);
    
    init_times(n_id) = data(1);
    odeint_times(n_id) = data(2);
end

fig = figure;
yyaxis left
hLine = plot(nums_particles, init_times, 'LineWidth', 2);
set(gca, 'FontSize', 20);
xlabel('$N$', 'Interpreter', 'latex');
set(gca, 'FontSize', 20);
ylabel('time of initialization', 'Interpreter', 'latex');
yyaxis right
hLine = plot(nums_particles, odeint_times, 'LineWidth', 2);
ylabel('time of $10^4$ integration steps', 'Interpreter', 'latex');

task_path = sprintf('%s(%d_%d_%0.2e)', ...
	task, ...
	start_observed_period, ...
	finish_observed_period, ...
	step);

new_dir = sprintf('%s/%s', figures_path, task_path);
mkdir(new_dir);
	
suffix = sprintf('np(var)_diss(%d_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)', ...
    diss_type, ...
    diss_gamma, ...
    E, ...
    U, ...
    J, ...
    drv_type, ...
    drv_ampl, ...
    drv_freq, ...
    drv_phase);

savefig(sprintf('%s/%s/run_times_from_np_%s.fig', figures_path, task_path, suffix));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/%s/run_times_from_np_%s.pdf', figures_path, task_path, suffix));