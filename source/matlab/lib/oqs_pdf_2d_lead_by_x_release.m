function pdf = oqs_pdf_2d_lead_by_x_release(pdf)

pdf.norm = zeros(pdf.x_num_points, 1);
for x_id = 1:pdf.x_num_points
    pdf.pdf(x_id, :) = pdf.pdf(x_id, :) / (pdf.inc_count(x_id) * pdf.y_bin_shift);
    pdf.norm(x_id) = sum(sum(pdf.pdf(x_id, :))) * pdf.y_bin_shift;
    pdf.pdf(x_id, :) = pdf.pdf(x_id, :) / max(pdf.pdf(x_id, :));
end

fprintf('max(pdf_norm) = %0.16e\n', max(pdf.norm));
fprintf('min(pdf_norm) = %0.16e\n', min(pdf.norm));

end
