function mtx = oqs_read_sparse_mtx(fn, N)

mtx_data = importdata(fn);
d_size = size(mtx_data, 1);

mtx = zeros(N,N);
for d_id = 1:d_size
    str = string(mtx_data(d_id));
    curr_data = sscanf(str, '%d\t%d\t(%e,%e)', 4);
    curr_row = curr_data(1) + 1;
    curr_col = curr_data(2) + 1;
    mtx(curr_row, curr_col) = curr_data(3) + 1i * curr_data(4);
end

end
