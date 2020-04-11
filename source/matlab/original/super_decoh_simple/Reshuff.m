function y = Reshuff(rho,N,M)

y = zeros(M,M);

for s = 1:M
    for sp = 1:M
        ii = mod(s,N);
        if (ii == 0)
            ii = N;
        end
        jj = 1 + (s - ii)/N;
        %
        kk = mod(sp,N);
        if (kk == 0)
            kk = N;
        end
        ll = 1 + (sp - kk)/N;
        %
        t=kk+N*(jj-1);
        tp=ll+N*(ii-1);
        y(s,sp)=rho(t,tp);
        
        
      %  t = ii + N*(kk - 1);
      %  tp = jj + N*(ll - 1);
       
    % y(s,sp) = rho(t,tp);
    end
end

end