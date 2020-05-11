function [alpha, coef, R2, yy, R2_, cnt] = oqs_powerlaw_regression(x, y, mn, mx)
    if nargin < 3
        mn = min(x);
    end
    if nargin < 4
        mx = max(x);
    end
    xx = x;
    y0 = y;
    ids = (y > 0) & (mn < x) & (x < mx);
    x = x(ids);
    y = y(ids);
    
    cnt = size(x);
    
    logx = log(x);
    logy = log(y);
    X = [ones(length(logx),1) logx];
    c = X \ logy;
    logyy = X * c;
    S1 = sum((logy - logyy).^2);
    S2 = sum((logy - mean(logy)).^2);
    R2 = 1 - S1 / S2;
    b = c(1);
    a = c(2);
    coef = exp(b);
    alpha = a;
    yy = coef * xx .^ alpha;
    R2_ = 1 - sum((y0 - yy).^2) / sum((y0 - mean(y0)).^2);
	
end