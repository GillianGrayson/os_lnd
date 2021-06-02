function y = Gen_Lindblad_from_Kolmogorov(seed, N)

imag1 = sqrt(-1.);

rng(seed);
MR = randn(N,'double');
MI = randn(N,'double');

MS = zeros(N,N);
MS  =MR + imag1 * MI;

MZ = zeros(N,N);
for j = 1:N
    for k = 1:N
        MZ(k,j) = MS(k,j)*conj(MS(k,j));
    end
end

M = N * N;


%*****************************************************************************************

X=zeros(M,M);
XX=zeros(M,M);
XX=1;

norm1=zeros(1,M);

for j = 1:N
    for k = 1:N
        s=N*(j-1)+k;
        X(s,s)=MZ(j,k);
        ttt=rand;
        norm1(s)=pi*(0.5-ttt);
    end
end


for j = 1:M
    for k = 1:M
        sss=norm1(j)-norm1(k);
        %XX(j,k)=sqrt(X(j,j)*X(k,k))*(cos(sss)+imag1*sin(sss));
    end
end

%****************************************************************************************
YR = randn(M,'double');
YI = randn(M,'double');
Y=YR+imag1*YI;
YY=Y*Y';

DD=zeros(M,M);

for j = 1:M
    DD(j,j)=1/sqrt(YY(j,j));
end

DELT1=YY*DD;
DELT=DD*DELT1;


%*************************************************************


rho_temp=times(DELT,XX);
y = sqrt(M)*rho_temp/trace(rho_temp);


