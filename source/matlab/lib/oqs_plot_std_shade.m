function fig = oqs_plot_std_shade(data)

amean = nanmean(data.data, 2);
astd = nanstd(data.data,[], 2); 
fig = figure;
h = fill([data.xs fliplr(data.xs)],[amean+astd fliplr(amean-astd)], 'FaceAlpha', data.alpha, 'linestyle', 'none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold all;
h = plot(data.xs, amean, 'LineWidth', 2);
legend(h, data.legend);
set(gca, 'FontSize', 30);
xlabel(data.x_label, 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel(data.y_label, 'Interpreter', 'latex');

end
