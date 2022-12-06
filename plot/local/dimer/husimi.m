function y = husimi(theta, phi, rho)

sizes = size(rho);
N = sizes(1)-1;
Ngrid = length(theta);
Mgrid = length(phi);
cos_theta = cos(theta/2);
sin_theta = sin(theta/2);
exp_iphi = exp(sqrt(-1)*phi);

y = zeros(Ngrid,Mgrid);

tmp = zeros(N);
for n=0:N
    tmp(n+1) = sqrt(nchoosek(N,n));
end

for i=1:Ngrid
    for j=1:Mgrid
        for n=0:N
            coherent(n+1,1) = tmp(n+1)*cos_theta(i)^n*(sin_theta(i)*exp_iphi(j))^(N-n);
        end
        y(i,j) = (coherent')*rho*coherent;
    end
end

y=y/sum(sum(y));

end