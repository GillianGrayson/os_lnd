function y = St_gen(num,M,Nr)

imag1 = sqrt(-1.);

%%% Generate X-matrix from Ginibre ensemble
seed1 = num;
rng(seed1, 'twister');
temp_y = randn(M,'double');

seed2 = Nr + 1 + num;
rng(seed2, 'twister');
temp_z = randn(M,'double');

X = (temp_y + imag1*temp_z)/sqrt(2);

%%% calculate rho-matrix
y = X*X';
y = sqrt(M)*y/trace(y);
end
