function pdf = oqs_pdf_1d_log_update(pdf, data)

for d_id = 1 : size(data, 1)
    
    x = data(d_id, 1);

    if ((x >= pdf.x_bin_s) && (x <= pdf.x_bin_f))
        x_id = floor((log10(x) - log10(pdf.x_bin_s)) * pdf.x_num_bins / (log10(pdf.x_bin_f) - log10(pdf.x_bin_s) + 1e-10)) + 1;
		pdf.pdf(x_id) = pdf.pdf(x_id) + 1;
        pdf.inc_count = pdf.inc_count + 1;
    else
        pdf.not_inc_count = pdf.not_inc_count + 1;
    end
end

end
