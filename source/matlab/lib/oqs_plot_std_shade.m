function fig = oqs_plot_std_shade(data)

amean = nanmean(data.data, 2);
astd = nanstd(data.data,[], 2);
fig = figure;
h = plot(data.xs, amean, 'LineWidth', 2);
legend(h, data.legend);
set(gca, 'FontSize', 30);
xlabel(data.x_label, 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel(data.y_label, 'Interpreter', 'latex');
hold all;
color = get(h,'Color');
h = fill(vertcat(data.xs, flip(data.xs)), vertcat(amean+astd, flip(amean-astd)), color, 'FaceAlpha', data.alpha, 'LineStyle', 'none');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';


end
