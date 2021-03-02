

thery_data_path = 'E:/YandexDisk/Work/os_lnd/data';
lind_3D = importdata(sprintf('%s/Kolmogorov_density_complex.dat', thery_data_path));
lind_real = importdata(sprintf('%s/Kolmogorov_density_real_eigenvalues.dat', thery_data_path));
lind_strip = importdata(sprintf('%s/Kolmogorov_density_complex_x_section.dat', thery_data_path));

figure;
h = plot(lind_real(:, 1), 0.97 * lind_real(:, 2), 'LineWidth', 2);
legend(h, 'real eigenvalues', 'interpeter', 'latex')
hold all;
h = plot(lind_strip(:, 1), 0.97 * 0.390 * sqrt(lind_strip(:, 2)), 'LineWidth', 2);
legend(h, 'y=0', 'interpeter', 'latex')
hold all;

pdf = pdf2d_quantum.pdf';

c_lim = 0.05;
width = 3;
for c_id = 1:size(pdf, 2)
    
    if (c_id == 7)
       ololo = 1; 
    end
    
    start_row = 1;
    for r_id = 1:101
        if pdf(r_id, c_id) > c_lim
            start_row = r_id;
            break;
        end
    end
    
    if start_row < 101
        mean_val = mean(pdf(start_row : 101 - width, c_id));
        replace = ones(2*width + 1, 1) * mean_val + randn(2*width + 1, 1) * 0.01 * mean_val;
        pdf(101 - width : 101 + width, c_id) = replace;
    end
end

pdf(101, :) = 0.5 * (pdf(100, :) + pdf(102, :));
ys = pdf(92, :);
xs = pdf2d_quantum.x_bin_centers;
h = plot(xs, 0.555 * sqrt(ys), 'LineWidth', 2);
legend(h, 'Sampling', 'interpeter', 'latex')
set(gca, 'FontSize', 30);
xlabel('$Re(\lambda)$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('PDF', 'Interpreter', 'latex');
