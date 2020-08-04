function pdf = oqs_pdf_2d_lead_by_x_update(pdf, data, x_id)

for y_id = 1:size(data, 1)
    y = data(y_id);
    if ((y >= pdf.y_bin_s) && (y <= pdf.y_bin_f))
        y_index = floor((y - pdf.y_bin_s) / (pdf.y_bin_shift + 1e-10)) + 1;
        pdf.pdf(x_id, y_index) = pdf.pdf(x_id, y_index) + 1;
        pdf.inc_count(x_id) = pdf.inc_count(x_id) + 1;
    else
        pdf.not_inc_count(x_id) = pdf.not_inc_count(x_id) + 1;
    end
end

end
