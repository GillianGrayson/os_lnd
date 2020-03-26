function pdf = oqs_pdf_2d_lead_by_x_setup(pdf)

pdf.y_bin_shift = (pdf.y_bin_f - pdf.y_bin_s) / pdf.y_num_bins;

pdf.y_bin_centers = linspace(...
    pdf.y_bin_s + 0.5 * pdf.y_bin_shift, ...
    pdf.y_bin_f - 0.5 * pdf.y_bin_shift, ...
    pdf.y_num_bins);

pdf.pdf = zeros(pdf.x_num_points, pdf.y_num_bins);
pdf.inc_count = zeros(pdf.x_num_points, 1);
pdf.not_inc_count =  zeros(pdf.x_num_points, 1);

end
