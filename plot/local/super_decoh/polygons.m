clear all;

path = 'C:/Users/user/Google Drive/os_lnd/data/superdecoh';

data = importdata(sprintf('%s/borderline_quantum.dat', path));
x1 = data(:, 1);
y1 = data(:, 2);
x2 = flip(x1);
y2 = -flip(y1);

x = vertcat(x1, x2);
y = vertcat(y1, y2);

pgon1 = polyshape(x, y);

%figure;
%plot(pgon1);

A = area(pgon1);

a = get(gca,'Children');
xdata = get(a, 'XData')';
ydata = get(a, 'YData')';
zdata = get(a, 'CData')';

xn1 = xdata;
yn1 = zeros(size(zdata, 2), 1);
for x_id = 1:size(zdata, 2)
    is_found = 0;
    for y_id = 1:size(zdata, 1)
        if zdata(x_id, y_id) > 0.001
            yn1(x_id) = ydata(y_id);
            is_found = 1;
            break;
        end
        if is_found == 0
            yn1(x_id) = 0;
        end
    end
end

yn1 = yn1 * -1.0;
xn2 = flip(xn1);
yn2 = -flip(yn1);

xn = vertcat(xn1, xn2);
yn = vertcat(yn1, yn2);
pgon2 = polyshape(xn, yn);

hold all;
plot(pgon2);

B = area(pgon2);


