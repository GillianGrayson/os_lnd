function y = Ggen(seed, M, N)

%%% Generate X-matrix
temp_y = randn(M,'double');
temp_z = randn(M,'double');

X = (temp_y + 1i*temp_z)/2.;

%%% calculate G-matrix
y = X*X';
y = N*y/trace(y);
end
