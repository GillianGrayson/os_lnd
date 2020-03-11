function fig = oqs_pdf_1d_log_plot(pdf)

fig = figure;
plot(pdf.x_bin_centers, pdf.pdf + 1e-10, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel(pdf.x_label, 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$PDF$', 'Interpreter', 'latex');
set(gca,'XScale','log');
set(gca,'YScale','log');
hold all;

end
