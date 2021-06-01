function y = Ginibre(seed, N)

rng(seed);

re = randn(N,'double');
im = randn(N,'double');

X = (re + 1i*im)/2.;

y = sqrt(N) * X / trace(X);

end
