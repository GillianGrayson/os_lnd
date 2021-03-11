clear all;
path = '../../../source/cpp/os_lnd/os_lnd';

serial_start = 1;
serial_shift = 1;
serial_num = 50;

Delta = 1.0;
W = 0.25;

mu = 0.001;

drv_type = 0;

T1 = 0.2;
T2 = 0.4;

quantity_index = 0;

num_times_per_seed = 3;
s_time_id = 3;
f_time_id = 3;

num_spins = [3,4,5,6,7]';
vars = zeros(size(num_spins, 1), 1);

opacity = 0.7;
legend_location = 'NorthEast';

for n_id = 1:size(num_spins, 1)
    
    ns = num_spins(n_id);
    
    suffix = sprintf('serial(%0.4f_%0.4f_%d)_ns(%d)_prm(%0.4f_%0.4f_%0.4f_%d_%0.4f_%0.4f)_j(%d)', ...
        serial_start, ...
        serial_shift, ...
        serial_num, ...
        ns, ...
        Delta, ...
        W, ...
        mu, ...
        drv_type, ...
        T1, ...
        T2, ...
        quantity_index);
    
    fn = sprintf('%s/serial_znd_%s.txt', path, suffix);
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
R2 = lm.Rsquared.Ordinary
RMSE = lm.RMSE;
intercept = lm.Coefficients{'(Intercept)', 'Estimate'};
slope = lm.Coefficients{'N', 'Estimate'};

fig = figure;
propertyeditor('on');
grid on;

h = scatter(num_spins, vars, 250, 'o', 'LineWidth',  1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerEdgeAlpha', opacity, 'MarkerFaceAlpha', opacity);
tmp = sprintf('$\\Delta=%0.2f$, $W=%0.2f$: $\\alpha=%0.3f$', Delta, W, slope);
legend(h, tmp, 'Interpreter', 'latex');
hold all;

x_shift = max(logN) - min(logN);
x_fit_log = linspace(min(logN) - 0.1 * x_shift,  max(logN) + 0.1 * x_shift, 1000);
y_fit_log = intercept + x_fit_log * slope;
x_fit = 10.^x_fit_log;
y_fit = 10.^y_fit_log;
h = plot(x_fit, y_fit, 'LineWidth', 2, 'Color', 'red');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold all;

set(gca, 'FontSize', 40);
xlabel('$N$', 'Interpreter', 'latex');
set(gca, 'FontSize', 40);
ylabel('$j/\mu$', 'Interpreter', 'latex');
legend(gca,'off');
legend('Location', legend_location, 'NumColumns', 1, 'Interpreter', 'latex');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
box on;


