function fig = oqs_pdf_2d_plot(pdf)

fig = figure;
imagesc(pdf.x_bin_centers, pdf.y_bin_centers, pdf.pdf');
set(gca, 'FontSize', 30);
xlabel(pdf.x_label, 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel(pdf.y_label, 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'PDF', 'FontSize', 33, 'interpreter','latex');
set(h, 'TickLabelInterpreter', 'latex');
set(gca,'YDir','normal');
ax = gca;
set(ax,'TickLabelInterpreter','Latex')
hold all;

end
