function pdf = oqs_pdf_1d_log_setup(pdf)

pdf.x_bin_borders = zeros(pdf.x_num_bins + 1, 1);
pdf.x_bin_centers = zeros(pdf.x_num_bins, 1);

pdf.x_bin_shift_log = (log10(pdf.x_bin_f) - log10(pdf.x_bin_s)) / pdf.x_num_bins;

for bin_id = 1 : pdf.x_num_bins + 1
    pdf.x_bin_borders(bin_id) = 10.^(log10(pdf.x_bin_s) + (bin_id - 1) * pdf.x_bin_shift_log);
    if (bin_id <= pdf.x_num_bins)
        pdf.x_bin_centers(bin_id) = 10.^(log10(pdf.x_bin_s) + (bin_id - 0.5) * pdf.x_bin_shift_log);
    end
end
pdf.x_bin_diff = diff(pdf.x_bin_borders);

pdf.pdf = zeros(pdf.x_num_bins, 1);
pdf.inc_count = 0;
pdf.not_inc_count = 0;

end
