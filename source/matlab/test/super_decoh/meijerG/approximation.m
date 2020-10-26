clear all;
x = linspace(-2, 2, 500);
z = x.^2 / 4;

y1 = meijerG([1, 2], [], [], [-0.5, 0.5],  z)';
aaa = meijerG(1, 2, -0.5, 0.5,  z)';
y2 = zeros(size(y1, 1), 1);
for y_id = 1:size(y1, 1)
    y2(y_id) = y1(y_id) * abs(x(y_id) / pi);
end


ololo = 1;

hold all;
plot(x, y2);