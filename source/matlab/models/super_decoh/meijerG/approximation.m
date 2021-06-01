clear all;
x = linspace(-2, 2, 500)';
x2 = x.^2;
tmp = 4.0 ./ x2;
arg = 1. - 4./(x.^2);

second = ellipticE(arg);
first = ellipke(arg);


y = 2 * abs(x) * ((4.+x.^2) * second - 8. * first)/(3 * pi.^2);

y1 = meijerG([1, 2], [], [], [-0.5, 0.5],  z)';
aaa = meijerG(1, 2, -0.5, 0.5,  z)';
y2 = zeros(size(y1, 1), 1);
for y_id = 1:size(y1, 1)
    y2(y_id) = y1(y_id) * abs(x(y_id) / pi);
end


ololo = 1;

hold all;
plot(x, y2);