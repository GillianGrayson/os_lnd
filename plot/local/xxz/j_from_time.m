clear all;
path = '../../../source/cpp/os_lnd/os_lnd';
fig_path = 'E:/YandexDisk/Work/os_lnd/figures/xxz/local';

fontSize = 30;
opacity = 1;
lineWidth = 1;
offLegend = 0;

num_spins = 7;

seed = 10;

Delta = 1.96;
W = 0.98;

mu = 0.001;

drv_type = 0;

ampl = 2.5;
freq = 1.0;
phase = 0.0;

j_id = 0;

num_periods = 51;
T = 2 * pi / freq;
final_time = T * num_periods;

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

fn = sprintf('%s/times_times(0.00e+00_%0.2e)_%s_seed(%d).txt', path, final_time, suffix, seed);
times = importdata(fn);
times = times / T;

fig_diag = figure;
propertyeditor('on');

for s_id = 1 : num_spins - 1
    fn = sprintf('%s/j_%d_times(0.00e+00_%0.2e)_%s_seed(%d).txt', path, s_id - 1, final_time, suffix, seed);
    var = importdata(fn);
    %var = log10(var);
    
    h = plot(times, var, 'LineWidth', lineWidth);
    h.Color = [h.Color opacity];
    if (offLegend == 0)
        legend(h, sprintf('k = %d', s_id - 1));
    else
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    set(gca, 'FontSize', fontSize);
    xlabel('$t/T$', 'Interpreter', 'latex');
    ylabel('$\log_{10}j_k$', 'Interpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex')
    box on;
    hold all;
    
end
set(gca, 'XScale', 'log')
legend(gca,'off');
legend('Location', 'SouthEast', 'NumColumns', 1, 'Interpreter', 'latex');
legend('FontSize', 20);
xlim([min(times), max(times)])

fn = sprintf('%s/j_from_time_%s_seed(%d)', fig_path, suffix, seed);
oqs_save_fig(gcf, fn);

