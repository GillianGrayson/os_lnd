function y = Hgen(seed, N)

%%% Generate GUE-matrix X
temp = randn(N,'double');

X = temp + 1i*temp;

%%% calculate G-matrix 
y = (X + X')/2.;
y2 = y*y;
y2_tr = trace(y2);
y = y/sqrt(y2_tr); 
end
