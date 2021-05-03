clear all;
addpath('../../../source/matlab/lib')

num_spins = [5, 6, 7]';
mu = 0.001;
drv_type = 0;
ampl = 1.0;
freq = 2 * pi;
phase = 0.0;

j_index = 0;

Delta = 1.0;
W = 0.5;

data_path = sprintf('/data/biophys/denysov/yusipov/os_lnd/serial/xxz/odeint');
figures_path = sprintf('/home/denysov/yusipov/os_lnd/figures/xxz');
if ~exist(figures_path, 'dir')
    mkdir(figures_path)
end

serial_start = 1;
serial_shift = 1;
serial_num = 100;

num_times_per_seed = 22;
s_time_id = 2;
f_time_id = 22;

vars = zeros(size(num_spins, 1), 1);
for n_id = 1:size(num_spins, 1)
    
    ns = num_spins(n_id);
    
    local_path = sprintf('n_%d/params_%0.4f_%0.4f_%0.4f_%d_%0.4f_%0.4f_%0.4f/seeds_%d_%d_%d', ...
        ns, ...
        Delta, ...
        W, ...
        mu, ...
        drv_type, ...
        ampl, ...
        freq, ...
        phase, ...
        serial_start, ...
        serial_shift, ...
        serial_num);
    
    suffix = sprintf('serial(%0.4f_%0.4f_%d)_ns(%d)_prm(%0.4f_%0.4f_%0.4f_%d_%0.4f_%0.4f_%0.4f)', ...
        serial_start, ...
        serial_shift, ...
        serial_num, ...
        ns, ...
        Delta, ...
        W, ...
        mu, ...
        drv_type, ...
        ampl, ...
        freq, ...
        phase);
    
    fn = sprintf('%s/%s/serial_j_%d_%s.txt', data_path, local_path, j_index, suffix);
    data = importdata(fn);
    
    curr_vars = zeros(serial_num, 1);
    for s_id = 1:serial_num
        start_id = (s_id - 1) * num_times_per_seed;
        tmp = data(start_id + s_time_id : start_id + f_time_id);
        curr_vars(s_id) = mean(tmp);
    end
    
    vars(n_id) = mean(curr_vars);
end
        
vars = vars / mu;
logVars = log10(vars);
logN = log10(num_spins);

T = table(logN, logVars, 'VariableNames', {'N', 'j'});
lm = fitlm(T, 'j~N');
R2 = lm.Rsquared.Ordinary;
RMSE = lm.RMSE;
intercept = lm.Coefficients{'(Intercept)', 'Estimate'};
slope = lm.Coefficients{'N', 'Estimate'};
log_x_shift = max(logN) - min(logN);
log_x_s =  min(logN) - 0.2 * log_x_shift;
log_x_f =  max(logN) + 0.2 * log_x_shift;
log_y_s = intercept + slope * log_x_s;
log_y_f = intercept + slope * log_x_f;
x_s = 10.^log_x_s;
x_f = 10.^log_x_f;
y_s = 10.^log_y_s;
y_f = 10.^log_y_f;

fig = figure;
h = plot(num_spins, vars, 'o');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold all;
h = plot([x_s, x_f], [y_s, y_f], 'LineWidth', 2);
legend(h, sprintf('$\\gamma$ = %0.2e', slope));
set(gca, 'FontSize', 30);
xlabel('$L$', 'Interpreter', 'latex');
ylabel('$j_k / \mu$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex')
box on;
hold all;
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
legend(gca,'off');
legend('Location', 'NorthEast', 'NumColumns', 1, 'Interpreter', 'latex');
legend('FontSize', 20);

suffix = sprintf('serial(%0.4f_%0.4f_%d)_prm(%0.4f_%0.4f_%0.4f_%d_%0.4f_%0.4f_%0.4f)_j(%d)', ...
    serial_start, ...
    serial_shift, ...
    serial_num, ...
    Delta, ...
    W, ...
    mu, ...
    drv_type, ...
    ampl, ...
    freq, ...
    phase, ...
    j_index);

fn_fig = sprintf('%s/x(N)_y(gamma)%s', figures_path, suffix);
oqs_save_fig(fig, fn_fig)
