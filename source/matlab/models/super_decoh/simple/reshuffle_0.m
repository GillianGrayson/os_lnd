function res = reshuffle_0(mtx, N)
N2 = N * N;
res = zeros(N2);
for i = 1:N
    for j = 1:N
        for k = 1:N
            for m = 1:N
                origin_row = (i - 1) * N + k;
                origin_col = (j - 1) * N + m;
                
                gamma_row = (m - 1) * N + k;
                gamma_col = (j - 1) * N + i;
                
                res(gamma_row, gamma_col) = mtx(origin_row, origin_col);
            end
        end
    end
end
end