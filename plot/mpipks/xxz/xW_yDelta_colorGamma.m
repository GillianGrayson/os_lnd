clear all;
addpath('../../../source/matlab/lib')

num_spins = [5, 6, 7]';
mu = 0.001;
drv_type = 0;
ampl = 2.5;
freq = 1.0;
phase = 0.0;

j_type = 'mean';
j_index = 0;

Deltas = [0.0:0.04:2.00]';
numDeltas = size(Deltas, 1);
Ws = [0.0:0.02:1.00]';
numWs = size(Ws, 1);

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

gammas = zeros(numWs, numDeltas);

js_all = {};
for n_id = 1:size(num_spins, 1)
    js_all{n_id} = zeros(numWs, numDeltas);
end

for W_id = 1:numWs
    W_id = W_id
    for D_id = 1:numDeltas
        Delta = Deltas(D_id);
        W = Ws(W_id);
        
        vars = zeros(size(num_spins, 1), 1);
        for n_id = 1:size(num_spins, 1)
    
            ns = num_spins(n_id);
            
            if strcmp(j_type, 'single')
                
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
                if isfile(fn)
                    data = importdata(fn);
                else
                    error('File does not exist: %s\n', fn);
                end
                
                
                curr_vars = zeros(serial_num, 1);
                for s_id = 1:serial_num
                    start_id = (s_id - 1) * num_times_per_seed;
                    tmp = data(start_id + s_time_id : start_id + f_time_id);
                    curr_vars(s_id) = mean(tmp);
                end
                
                vars(n_id) = mean(curr_vars);
                
            elseif strcmp(j_type, 'mean')
                
                curr_js_means = zeros(ns-1, 1);
                for curr_j = 1:ns-1
                
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
                    
                    fn = sprintf('%s/%s/serial_j_%d_%s.txt', data_path, local_path, curr_j-1, suffix);
                    if isfile(fn)
                        data = importdata(fn);
                    else
                        error('File does not exist: %s\n', fn);
                    end

                    curr_vars = zeros(serial_num, 1);
                    for s_id = 1:serial_num
                        start_id = (s_id - 1) * num_times_per_seed;
                        tmp = data(start_id + s_time_id : start_id + f_time_id);
                        curr_vars(s_id) = mean(tmp);
                    end
                    
                    curr_js_means(curr_j) = mean(curr_vars);
                end
                
                vars(n_id) = mean(curr_js_means);
                
            end
        end
        
        vars = vars / mu;
        logVars = log10(vars);
        
        for n_id = 1:size(num_spins, 1)
            js_all{n_id}(W_id, D_id) = vars(n_id);
        end
        
        logN = log10(num_spins);
        
        T = table(logN, logVars, 'VariableNames', {'N', 'j'});
        lm = fitlm(T, 'j~N');
        R2 = lm.Rsquared.Ordinary;
        RMSE = lm.RMSE;
        intercept = lm.Coefficients{'(Intercept)', 'Estimate'};
        slope = lm.Coefficients{'N', 'Estimate'};
        
        gammas(W_id, D_id) = -slope;
    end
end


if strcmp(j_type, 'single')
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
elseif strcmp(j_type, 'mean')
    suffix = sprintf('serial(%0.4f_%0.4f_%d)_prm(%0.4f_%0.4f_%0.4f_%d_%0.4f_%0.4f_%0.4f)_j(mean)', ...
    serial_start, ...
    serial_shift, ...
    serial_num, ...
    Delta, ...
    W, ...
    mu, ...
    drv_type, ...
    ampl, ...
    freq, ...
    phase);
end

for n_id = 1:size(num_spins, 1)
    fig = figure;
    imagesc(Ws, Deltas, js_all{n_id}');
    set(gca, 'FontSize', 30);
    xlabel('$W$', 'Interpreter', 'latex');
    set(gca, 'FontSize', 30);
    ylabel('$\Delta$', 'Interpreter', 'latex');
    colormap jet;
    h = colorbar;
    set(gca, 'FontSize', 30);
    title(h, '$j/\mu$', 'FontSize', 33, 'interpreter','latex');
    set(gca,'YDir','normal');
    hold all;
    fn_fig = sprintf('%s/x(W)_y(Delta)_color(j)_n(%d)_%s', figures_path, num_spins(n_id), suffix);
    oqs_save_fig(fig, fn_fig)
end

fig = figure;
imagesc(Ws, Deltas, gammas');
set(gca, 'FontSize', 30);
xlabel('$W$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\Delta$', 'Interpreter', 'latex');
colormap jet;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$\gamma$', 'FontSize', 33, 'interpreter','latex');
set(gca,'YDir','normal');
hold all;
fn_fig = sprintf('%s/x(W)_y(Delta)_color(gamma)_%s', figures_path, suffix);
oqs_save_fig(fig, fn_fig)
