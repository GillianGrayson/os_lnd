function fig = oqs_pdf_1d_plot(pdf)

fig = figure;
plot(pdf.x_bin_centers, pdf.pdf, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel(pdf.x_label, 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$PDF$', 'Interpreter', 'latex');
hold all;

end
