function pdf = oqs_pdf_2d_lead_by_x_release(pdf, is_smooth, is_normalize)

pdf.norm = zeros(pdf.x_num_points, 1);
for x_id = 1:pdf.x_num_points
    pdf.pdf(x_id, :) = pdf.pdf(x_id, :) / (pdf.inc_count(x_id) * pdf.y_bin_shift);
	if is_smooth > 0
		pdf.pdf(x_id, :) = smoothdata(pdf.pdf(x_id, :));
	end
	if is_normalize > 0
		pdf.pdf(x_id, :) = pdf.pdf(x_id, :) / max(pdf.pdf(x_id, :));
	end
	pdf.norm(x_id) = sum(sum(pdf.pdf(x_id, :))) * pdf.y_bin_shift;
end

fprintf('max(pdf_norm) = %0.16e\n', max(pdf.norm));
fprintf('min(pdf_norm) = %0.16e\n', min(pdf.norm));

end
