function [i_max, j_max] = oqs_find_range(x, y)
    i_max = 1; j_max = 1;
    ids = y > 0;
    x = x(ids);
    y = y(ids);
    len = 1.0;
    R2_max = 0;
    for i = 1:size(x)
        for j = i:size(x)
            [~, ~, R2_pow, ~, ~, cnt] = powerlaw_regression(x, y, x(i), x(j));
            if cnt < 10
                continue;
            end
            curlen = log10(x(j)) - log10(x(i));
            if R2_pow > 0.98 &&  curlen > 1.0
                 curlen = curlen + 120 * (R2_pow - 0.98);
                if abs(len - curlen) <= 0.01
                    if R2_pow > R2_max
                        len = curlen;
                        R2_max = R2_pow;
                        i_max = x(i);
                        j_max = x(j);
                    end
                elseif len + 1 < curlen
                    len = curlen;
                    R2_max = R2_pow;
                    i_max = x(i);
                    j_max = x(j);
                end
            end
        end
    end 
end

