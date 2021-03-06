function fig = oqs_pdf_2d_lead_by_x_plot(pdf)

fig = figure;
imagesc(pdf.xs, pdf.y_bin_centers, pdf.pdf');
set(gca, 'FontSize', 30);
xlabel(pdf.x_label, 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel(pdf.y_label, 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$PDF$', 'FontSize', 33, 'interpreter','latex');
set(gca,'YDir','normal');
hold all;

end
