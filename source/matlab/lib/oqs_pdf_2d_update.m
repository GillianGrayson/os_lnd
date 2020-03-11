function pdf = oqs_pdf_2d_update(pdf, data)

for d_id = 1 : size(data, 1)
    
    x = data(d_id, 1);
    y = data(d_id, 2);
    
    if ((x >= pdf.x_bin_s) && (x <= pdf.x_bin_f) && (y >= pdf.y_bin_s) && (y <= pdf.y_bin_f))
        x_id = floor((x - pdf.x_bin_s) / (pdf.x_bin_shift + 1e-10)) + 1;
        y_id = floor((y - pdf.y_bin_s) / (pdf.y_bin_shift + 1e-10)) + 1;
        
        pdf.pdf(x_id, y_id) = pdf.pdf(x_id, y_id) + 1;
        pdf.inc_count = pdf.inc_count + 1;
    else
        pdf.not_inc_count = pdf.not_inc_count + 1;
    end
end

end
