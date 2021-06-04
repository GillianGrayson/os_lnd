clear all;
addpath('../../../source/matlab/lib')

task = 'eigen_dense';

mem_size = 2;
non_zeros_size = 1;
run_times_size = 4;

path = sprintf('/data/biophys/denysov/yusipov/os_lnd/regular/lfk/%s', task);

Ns = [100]';
total_num_evals = 100000;
num_seeds = zeros(size(Ns, 1), 1);

for N_id = 1:size(Ns, 1)
    
    N = Ns(N_id);
    N2 = N * N;
    
    num_seeds(N_id) = 1000;
    
    fprintf('N = %d\n', N);
    fprintf('num_seeds = %d\n', num_seeds(N_id));
    
    
    mem_infos_all = zeros(mem_size * num_seeds(N_id), 1);
    non_zeros_all = zeros(non_zeros_size * num_seeds(N_id), 1);
    run_times_all = zeros(run_times_size * num_seeds(N_id), 1);
    rho_evals_all = zeros(N * num_seeds(N_id), 1);
    lind_evals_all = zeros(N2 * num_seeds(N_id), 1);
    
    for seed_id = 1:num_seeds(N_id)
        
        seed = seed_id
        
        suffix = sprintf('N(%d)_seed(%d)', ...
            N, ...
            seed);
        
        prefix = sprintf('%s/N_%d/seed_%d', ...
            path, ...
            N, ...
            seed);
        
        fn = sprintf('%s/mem_info_%s.txt', ...
            prefix, ...
            suffix);
        mem_info_data = importdata(fn);
        mem_infos_all((seed_id - 1) * mem_size + 1 : seed_id * mem_size) = max(mem_info_data);
        
        fn = sprintf('%s/non_zeros_parts_%s.txt', ...
            prefix, ...
            suffix);
        non_zero_data = importdata(fn);
        non_zeros_all((seed_id - 1) * non_zeros_size + 1 : seed_id * non_zeros_size) = non_zero_data;
        
        fn = sprintf('%s/run_times_%s.txt', ...
            prefix, ...
            suffix);
        run_times_data = importdata(fn);
        run_times_all((seed_id - 1) * run_times_size + 1 : seed_id * run_times_size) = max(run_times_data);
        
        fn = sprintf('%s/rho_evals_%s.txt', ...
            prefix, ...
            suffix);
        rho_evals_data = importdata(fn);
        evals = zeros(N, 1);
        for str_id = 1:N
            str = string(rho_evals_data(str_id));
            data2d = sscanf(str, '(%e,%e)', 2);
            evals(str_id) = data2d(1) + 1i * data2d(2);
        end
        rho_evals_all((seed_id - 1) * N + 1 : seed_id * N) = evals;
        
        fn = sprintf('%s/lindbladian_evals_%s.txt', ...
            prefix, ...
            suffix);
        lindbladian_evals_data = importdata(fn);
        evals = zeros(N2, 1);
        for str_id = 1:N2
            str = string(lindbladian_evals_data(str_id));
            data2d = sscanf(str, '(%e,%e)', 2);
            evals(str_id) = data2d(1) + 1i * data2d(2);
        end
        lind_evals_all((seed_id - 1) * N2 + 1 : seed_id * N2) = evals;
    end
    
    prefix = sprintf('N_%d', ...
        N);
    
    aggr_path = sprintf('%s/aggregator/%s', ...
        path, ...
        prefix);
    mkdir(aggr_path);
    
    suffix = sprintf('N(%d)_seeds(%d_%d_%d)', ...
        N, ...
        1, ...
        1, ...
        num_seeds(N_id));
    
    fn_txt = sprintf('%s/mem_info_%s.txt', aggr_path, suffix);
    fid = fopen(fn_txt,'wt');
    for x_id = 1:size(mem_infos_all, 1)
        fprintf(fid,'%0.16e\n', mem_infos_all(x_id));
    end
    fclose(fid);
    
    fn_txt = sprintf('%s/non_zeros_parts_%s.txt', aggr_path, suffix);
    fid = fopen(fn_txt,'wt');
    for x_id = 1:size(non_zeros_all, 1)
        fprintf(fid,'%0.16e\n', non_zeros_all(x_id));
    end
    fclose(fid);
    
    fn_txt = sprintf('%s/run_times_%s.txt', aggr_path, suffix);
    fid = fopen(fn_txt,'wt');
    for x_id = 1:size(run_times_all, 1)
        fprintf(fid,'%0.16e\n', run_times_all(x_id));
    end
    fclose(fid);
    
    fn_txt = sprintf('%s/rho_evals_%s.txt', aggr_path, suffix);
    fid = fopen(fn_txt,'wt');
    for x_id = 1:size(rho_evals_all, 1)
        fprintf(fid,'%0.16e %0.16e\n', real(rho_evals_all(x_id)), imag(rho_evals_all(x_id)));
    end
    fclose(fid);
    
    fn_txt = sprintf('%s/lindbladian_evals_%s.txt', aggr_path, suffix);
    fid = fopen(fn_txt,'wt');
    for x_id = 1:size(lind_evals_all, 1)
        fprintf(fid,'%0.16e %0.16e\n', real(lind_evals_all(x_id)), imag(lind_evals_all(x_id)));
    end
    fclose(fid);
    
end
