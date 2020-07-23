clear all;
addpath('../../../source/matlab/lib')

task = 'eigen_dense';

path = sprintf('/data/condmat/ivanchen/yusipov/os_lnd/mbl/%s', task);

ns = 8;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

Ws = [1, 20]';
U = 1.0;
J = 1.0;

seeds = linspace(1, 100, 100)';

num_Ws = size(Ws, 1);
num_seeds = size(seeds, 1);

N = (nchoosek(ns, ns/2));
N2 = N^2;

for W_id = 1:num_Ws
    
    W = Ws(W_id);
    fprintf('W = %0.16e\n', W);
    
    ees_all = zeros(size(seeds, 1), 1);
    imbalances_all = zeros(size(seeds, 1), 1);
    ratios_all = zeros(size(seeds, 1), 1);
    
    mem_infos_all = zeros(2 * size(seeds, 1), 1);
    non_zeros_all = zeros(2 * size(seeds, 1), 1);
    run_times_all = zeros(4 * size(seeds, 1), 1);
    
    energies_all = zeros(ns * size(seeds, 1), 1);
    
    rho_evals_all = zeros(N * size(seeds, 1), 1);
    
    lind_evals_all = zeros(N2 * size(seeds, 1), 1);
    
    for seed_id = 1:num_seeds
	
        seed = seeds(seed_id);
		fprintf('seed = %d\n', seed);
		
        suffix = sprintf('ns(%d)_seed(%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)', ...
			ns, ...
            seed, ...
            diss_type, ...
            diss_phase, ...
            diss_gamma, ...
            W, ...
            U, ...
            J);
        
        prefix = sprintf('%s/ns_%d/diss_%d_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/seed_%d', ...
            path, ...
            ns, ...
            diss_type, ...
            diss_phase, ...
            diss_gamma, ...
            W, ...
            U, ...
            J, ...
            seed);
        
        fn = sprintf('%s/ee_%s.txt', ...
            prefix, ...
            suffix);
        ee_data = importdata(fn);
        ees_all(seed_id) = ee_data;
     
        fn = sprintf('%s/imbalance_%s.txt', ...
            prefix, ...
            suffix);
        imbalance_data = importdata(fn);
        imbalances_all(seed_id) = imbalance_data;
        
        fn = sprintf('%s/ratio_%s.txt', ...
            prefix, ...
            suffix);
        ratio_data = importdata(fn);
        ratios_all(seed_id) = ratio_data;
        
        fn = sprintf('%s/mem_info_%s.txt', ...
            prefix, ...
            suffix);
        mem_info_data = importdata(fn);
        mem_infos_all((seed_id - 1) * 2 + 1 : seed_id * 2) = mem_info_data;
        
        fn = sprintf('%s/non_zeros_parts_%s.txt', ...
            prefix, ...
            suffix);
        non_zero_data = importdata(fn);
        non_zeros_all((seed_id - 1) * 2 + 1 : seed_id * 2) = non_zero_data;
        
        fn = sprintf('%s/run_times_%s.txt', ...
            prefix, ...
            suffix);
        run_times_data = importdata(fn);
        run_times_all((seed_id - 1) * 4 + 1 : seed_id * 4) = run_times_data;
        
        fn = sprintf('%s/energies_%s.txt', ...
            prefix, ...
            suffix);
        energies_data = importdata(fn);
        energies_all((seed_id - 1) * ns + 1 : seed_id * ns) = energies_data;
        
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
    
    prefix = sprintf('ns_%d/diss_%d_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f', ...
        ns, ...
        diss_type, ...
        diss_phase, ...
        diss_gamma, ...
        W, ...
        U, ...
        J);
    
    aggr_path = sprintf('%s/aggregator/%s', ...
        path, ...
        prefix);
    mkdir(aggr_path);
        
    suffix = sprintf('ns(%d)_seeds(%d_%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)', ...
        ns, ...
        1, ...
        1, ...
        num_seeds, ...
        diss_type, ...
        diss_phase, ...
        diss_gamma, ...
        W, ...
        U, ...
        J);
    
    fn_txt = sprintf('%s/ee_%s.txt', aggr_path, suffix);
    fid = fopen(fn_txt,'wt');
    for x_id = 1:size(ees_all, 1)
        fprintf(fid,'%0.16e\n', ees_all(x_id));
    end
    fclose(fid);
    
    fn_txt = sprintf('%s/imbalance_%s.txt', aggr_path, suffix);
    fid = fopen(fn_txt,'wt');
    for x_id = 1:size(imbalances_all, 1)
        fprintf(fid,'%0.16e\n', imbalances_all(x_id));
    end
    fclose(fid);
    
    fn_txt = sprintf('%s/ratio_%s.txt', aggr_path, suffix);
    fid = fopen(fn_txt,'wt');
    for x_id = 1:size(ratios_all, 1)
        fprintf(fid,'%0.16e\n', ratios_all(x_id));
    end
    fclose(fid);
    
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
    
    fn_txt = sprintf('%s/energies_%s.txt', aggr_path, suffix);
    fid = fopen(fn_txt,'wt');
    for x_id = 1:size(energies_all, 1)
        fprintf(fid,'%0.16e\n', energies_all(x_id));
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

