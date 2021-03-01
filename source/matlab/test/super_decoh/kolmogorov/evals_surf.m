f = figure;
pdf = pdf2d.pdf';

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

hLine = plot(xt, yt, 'LineWidth', 3, 'Color', 'cyan');


