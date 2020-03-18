function y = Ggen_alter(num,N,Nr)

imag1 = sqrt(-1.);

%%% Generate X-matrix
seed1 = num;
rng(seed1, 'twister');
temp_y = randn(N,'double');

seed2 = Nr + 1 + num;
rng(seed2, 'twister');
temp_z = randn(N,'double');

X = (temp_y + imag1*temp_z)/sqrt(2);

%%% calculate G-matrix
K = X*X';
TR = trace(K);
y = X*sqrt(N/TR);
end
