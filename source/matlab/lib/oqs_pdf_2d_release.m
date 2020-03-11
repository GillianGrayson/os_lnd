function pdf = oqs_pdf_2d_release(pdf)

pdf.pdf = pdf.pdf / (pdf.inc_count * pdf.x_bin_shift * pdf.y_bin_shift);
pdf.norm = sum(sum(pdf.pdf)) * pdf.x_bin_shift * pdf.y_bin_shift;
fprintf('pdf_norm = %0.16e\n', pdf.norm);

end
