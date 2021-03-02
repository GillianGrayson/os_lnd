f = figure;
pdf = pdf2d.pdf';

c_lim = 0.1;
width = 6;
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

pdf = smooth2a(pdf, 5, 5);

s = surf(pdf2d.x_bin_centers, pdf2d.y_bin_centers, pdf, 'FaceAlpha', 0.8);
set(gca, 'FontSize', 30);
xlabel(pdf2d.x_label, 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel(pdf2d.y_label, 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$PDF$', 'FontSize', 33, 'interpreter','latex');
set(gca,'YDir','normal');
hold all;
propertyeditor('on')

hold all;
theory_data_quantum = importdata(sprintf('%s/borderline_classical.dat', thery_data_path));
xt1 = theory_data_quantum(:, 1);
yt1 = theory_data_quantum(:, 2);
xt2 = flip(xt1);
yt2 = -flip(yt1);
xt = vertcat(xt1, xt2);
yt = vertcat(yt1, yt2);
hLine = plot3(xt, yt, zeros(size(xt, 1), 1), 'LineWidth', 3, 'Color', 'k');

hold all
slice_id = 101;
xs = pdf2d.x_bin_centers(slice_id) * ones(size(pdf2d.y_bin_centers, 2), 1);
ys = pdf2d.y_bin_centers;
zs = pdf(:, slice_id);
plot3(xs, ys, zs)

slice_id = 76;
xs = pdf2d.x_bin_centers(slice_id) * ones(size(pdf2d.y_bin_centers, 2), 1);
ys = pdf2d.y_bin_centers;
zs = pdf(:, slice_id);
plot3(xs, ys, zs)

slice_id = 126;
xs = pdf2d.x_bin_centers(slice_id) * ones(size(pdf2d.y_bin_centers, 2), 1);
ys = pdf2d.y_bin_centers;
zs = pdf(:, slice_id);
plot3(xs, ys, zs)

slice_id = 81;
xs = pdf2d.x_bin_centers;
ys = pdf2d.y_bin_centers(slice_id) * ones(size(pdf2d.x_bin_centers, 2), 1);
zs = pdf(slice_id, :);
plot3(xs, ys, zs)
zs_zeros = zeros(size(xs, 2), 1);
plot3(xs, ys, zs_zeros)

