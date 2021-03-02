

thery_data_path = 'E:/YandexDisk/Work/os_lnd/data';
lind_3D = importdata(sprintf('%s/Kolmogorov_density_complex.dat', thery_data_path));
lind_real = importdata(sprintf('%s/Kolmogorov_density_real_eigenvalues.dat', thery_data_path));
lind_strip = importdata(sprintf('%s/Kolmogorov_density_complex_x_section.dat', thery_data_path));

shift_y = 201;

xs = lind_3D(1:shift_y:end, 1);
ys = lind_3D(1:shift_y, 2);

pdf = zeros(size(xs, 1), size(ys, 1));
for x_id = 1:size(xs, 1)
    for y_id = 1:size(ys, 1)
        index = (x_id - 1) * size(ys, 1) + y_id;
        pdf(x_id, y_id) = lind_3D(index, 3);
    end
end

xs = xs(500:900);
pdf = pdf(500:900, :);

f = figure;
pdf = pdf';
s = surf(xs, ys, pdf, 'FaceAlpha', 0.8);
set(gca, 'FontSize', 30);
xlabel('$Re(\lambda)$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$Im(\lambda)$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$PDF$', 'FontSize', 33, 'interpreter','latex');
set(gca,'YDir','normal');
hold all;
propertyeditor('on')

slice_id = 600;
xs_curr = xs(slice_id) * ones(size(ys, 1), 1);
ys_curr = ys;
zs_curr = pdf(:, slice_id);
plot3(xs_curr, ys_curr, zs_curr)

slice_id = 700;
xs_curr = xs(slice_id) * ones(size(ys, 1), 1);
ys_curr = ys;
zs_curr = pdf(:, slice_id);
plot3(xs_curr, ys_curr, zs_curr)

slice_id = 800;
xs_curr = xs(slice_id) * ones(size(ys, 1), 1);
ys_curr = ys;
zs_curr = pdf(:, slice_id);
plot3(xs_curr, ys_curr, zs_curr)

slice_id = 81;
xs_curr = xs;
ys_curr = ys(slice_id) * ones(size(xs, 1), 1);
zs_curr = pdf(slice_id, :);
plot3(xs_curr, ys_curr, zs_curr)

hold all;
theory_data_quantum = importdata(sprintf('%s/borderline_classical.dat', thery_data_path));
xt1 = theory_data_quantum(:, 1);
yt1 = theory_data_quantum(:, 2);
xt2 = flip(xt1);
yt2 = -flip(yt1);
xt = vertcat(xt1, xt2);
yt = vertcat(yt1, yt2);
hLine = plot3(xt, yt, zeros(size(xt, 1), 1), 'LineWidth', 3, 'Color', 'k');
