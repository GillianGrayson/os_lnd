function y = Ggen(num,M,N)

imag1 = sqrt(-1.);

%%% Generate X-matrix
seed = num;
randn('state', seed);
temp_y = randn(M,'double');

seed = 1001 + num;
randn('state', seed);
temp_z = randn(M,'double');

X = (temp_y + imag1*temp_z)/2.;

%%% calculate G-matrix
y = X*X';
y = N*y/trace(y);
end
