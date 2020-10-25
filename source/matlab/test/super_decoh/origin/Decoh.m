function y = Decoh(rho_temp,M,p)

d = diag(rho_temp);
y = p*rho_temp;

for ii = 1:M
    y(ii,ii) = d(ii,1);
end

end