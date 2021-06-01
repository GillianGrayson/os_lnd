% ???? ?????? ???? ?????????? ???????????, ? ????? ???

imag_parts = [0.0]';

imag_lim = 1e-8;

for imag_id = 1:size(imag_parts, 1)
    imag_part = imag_parts(imag_id);
    
    passed_evals = zeros(size(size(all_evals, 2) * size(all_evals, 1), 1), 1);
    num_passed_values = 0;
    for seed = 1:size(all_evals, 2)
        curr_evals = all_evals(2:end, seed);
        curr_evals = sqrt(N) * (real(curr_evals) + 1) + 1i * sqrt(N) * imag(curr_evals);
        for eval_id = 1 : size(curr_evals, 1)
            if (imag(curr_evals(eval_id)) >= (imag_part - imag_lim)) && (imag(curr_evals(eval_id)) < (imag_part + imag_lim))
                num_passed_values = num_passed_values + 1;
                passed_evals(num_passed_values, 1) = real(curr_evals(eval_id));
            end
        end
    end
    
    passed_evals = passed_evals(1:num_passed_values);
    
    real_evals.x_num_bins = 201;
    real_evals.x_label = '$Re(\lambda)$';
    real_evals.x_bin_s = min(passed_evals) - 1e-16;
    real_evals.x_bin_f = max(passed_evals) + 1e-16;
    real_evals = oqs_pdf_1d_setup(real_evals);
    real_evals = oqs_pdf_1d_update(real_evals, passed_evals);
    real_evals = oqs_pdf_1d_release(real_evals);
    fig = oqs_pdf_1d_plot(real_evals);
end