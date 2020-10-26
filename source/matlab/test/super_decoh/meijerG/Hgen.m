function y = Hgen(seed, N)

rng(seed);
tmp = randn(N, 'double');
x = tmp + 1i * tmp;
y = (x + x') / sqrt(2);
y2 = y * y;
y2_tr = trace(y2);
y = sqrt(N)/2 * y/sqrt(y2_tr); 

end
