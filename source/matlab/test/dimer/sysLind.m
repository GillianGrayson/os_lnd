function dy = sysLind(t,x)
global N A0 w A phi0 Ps Pd

T=2*pi/w;

ft=sin(w*t+phi0);%t/T-floor(t/T);
%ft=2*(ft>0)-1;

dy=(Ps+Pd*A0*ft)*x;

end

