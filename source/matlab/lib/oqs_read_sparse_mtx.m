function mtx = oqs_read_sparse_mtx(fn, N)

mtx_data = importdata(fn);
d_size = size(mtx_data, 1);

mtx = zeros(N,N);
for d_id = 1:d_size
    curr_row = mtx_data(d_id, 1);
    curr_col = mtx_data(d_id, 2);
    mtx(curr_row, curr_col) = mtx_data(d_id, 3) + 1i * mtx_data(d_id, 4);
end

end
