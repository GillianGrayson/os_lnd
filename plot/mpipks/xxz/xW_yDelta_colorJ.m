clear all;
addpath('../../../source/matlab/lib')

num_spins = [6, 7]';
mu = 0.001;
drv_type = 0;
ampl = 0.0;
freq = 10 * pi;
phase = 0.0;
quantity_index = 0;

Deltas = linspace(0.0, 2.0, 51)';
numDeltas = size(Deltas, 1);
Ws = linspace(0.0, 1.0, 51)';
numWs = size(Ws, 1);

data_path = sprintf('/data/biophys/denysov/yusipov/os_lnd/serial/xxz/odeint');
figures_path = sprintf('/home/denysov/yusipov/os_lnd/figures/xxz');
if ~exist(figures_path, 'dir')
    mkdir(figures_path)
end

serial_start = 1;
serial_shift = 1;
serial_num = 100;

num_times_per_seed = 3;
s_time_id = 3;
f_time_id = 3;

js = zeros(numWs, numDeltas);

for W_id = 1:numWs
    W_id = W_id
    for D_id = 1:numDeltas
        Delta = Deltas(D_id);
        W = Ws(W_id);
        
        vars = zeros(size(num_spins, 1), 1);
        for n_id = 1:size(num_spins, 1)
    
            ns = num_spins(n_id);
            
            local_path = sprintf('n_%d/params_%0.4f_%0.4f_%0.4f_%d_%0.4f_%0.4f_%0.4f_%d/seeds_%d_%d_%d', ...
                ns, ...
                Delta, ...
                W, ...
                mu, ...
                drv_type, ...
                ampl, ...
                freq, ...
                phase, ...
                quantity_index, ...
                serial_start, ...
                serial_shift, ...
                serial_num);
            
            suffix = sprintf('serial(%0.4f_%0.4f_%d)_ns(%d)_prm(%0.4f_%0.4f_%0.4f_%d_%0.4f_%0.4f_%0.4f)_j(%d)', ...
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
                phase, ...
                quantity_index);
            
            fn = sprintf('%s/%s/serial_znd_%s.txt', data_path, local_path, suffix);
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
        
        js(W_id, D_id) = -slope;
    end
end

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
    quantity_index);

fig = figure;
imagesc(Ws, Deltas, js');
set(gca, 'FontSize', 30);
xlabel('$W$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\Delta$', 'Interpreter', 'latex');
colormap jet;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$j / \mu$', 'FontSize', 33, 'interpreter','latex');
set(gca,'YDir','normal');
hold all;
fn_fig = sprintf('%s/x(W)_y(Delta)_color(j)_%s', figures_path, suffix);
oqs_save_fig(fig, fn_fig)
