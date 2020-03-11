function pdf = oqs_pdf_2d_setup(pdf)

pdf.x_bin_shift = (pdf.x_bin_f - pdf.x_bin_s) / pdf.x_num_bins;
pdf.y_bin_shift = (pdf.y_bin_f - pdf.y_bin_s) / pdf.y_num_bins;

pdf.x_bin_centers = linspace(...
    pdf.x_bin_s + 0.5 * pdf.x_bin_shift, ...
    pdf.x_bin_f - 0.5 * pdf.x_bin_shift, ...
    pdf.x_num_bins);

pdf.y_bin_centers = linspace(...
    pdf.y_bin_s + 0.5 * pdf.y_bin_shift, ...
    pdf.y_bin_f - 0.5 * pdf.y_bin_shift, ...
    pdf.y_num_bins);

pdf.pdf = zeros(pdf.x_num_bins, pdf.y_num_bins);
pdf.inc_count = 0;
pdf.not_inc_count = 0;

end
