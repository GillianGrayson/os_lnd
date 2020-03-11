function pdf = oqs_pdf_1d_log_release(pdf)

pdf.norm = 0.0;
for bin_id = 1 : pdf.x_num_bins
    pdf.pdf(bin_id) = pdf.pdf(bin_id) / (pdf.inc_count * pdf.x_bin_diff(bin_id));
    pdf.norm = pdf.norm + pdf.pdf(bin_id) * pdf.x_bin_diff(bin_id);
end
fprintf('pdf_norm = %0.16e\n', pdf.norm);

end