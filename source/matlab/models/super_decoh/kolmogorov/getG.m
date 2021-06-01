function y = getG(N, seed)

rng(seed);
temp_y = randn(N,'double');
temp_z = randn(N,'double');
X = (temp_y + 1i * temp_z) / sqrt(2);
K = X*X';
TR = trace(K);
y = X*sqrt(N/TR);

end
