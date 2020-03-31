clear all;
addpath('../../../source/matlab/lib')

path = '/data/condmat/ivanchen/yusipov/os_lnd/mbl';
figures_path = '/home/ivanchen/yusipov/os_lnd/figures/mbl';

ns = 8;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

Ws = linspace(1, 20, 20)';
U = 1.0;
J = 1.0;

seeds = linspace(1, 100, 100)';

num_Ws = size(Ws, 1);
num_seeds = size(seeds, 1);

x_num_bins = 301;
y_num_bins = 301;
x_label = '$Re(\lambda)$';
y_label = '$Im(\lambda)$';

N2 = (nchoosek(ns, ns/2)).^2

for W_id = 1:num_Ws
    
    W = Ws(W_id);
    fprintf('W = %0.16e\n', W);
    
    pdf2d.x_num_bins = x_num_bins;
    pdf2d.y_num_bins = y_num_bins;
    pdf2d.x_label = x_label;
    pdf2d.y_label = y_label;
    
    all_evals = zeros((N2 - 1) * size(seeds, 1), 1);
    
    for seed_id = 1:num_seeds
	
        seed = seeds(seed_id);
		
        suffix = sprintf('ns(%d)_seed(%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)', ...
			ns, ...
            seed, ...
            diss_type, ...
            diss_phase, ...
            diss_gamma, ...
            W, ...
            U, ...
            J);
        
        fn = sprintf('%s/eigen_dense/ns_%d/diss_%d_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/seed_%d/lindbladian_evals_%s.txt', ...
            path, ...
            ns, ...
            diss_type, ...
            diss_phase, ...
            diss_gamma, ...
            W, ...
            U, ...
            J, ...
            seed, ...
            suffix);
        
        evals = zeros(N2, 1);
        data = importdata(fn);
        for str_id = 1:N2
            str = string(data(str_id));
            data2d = sscanf(str, '(%e,%e)', 2);
            evals(str_id) = data2d(1) + 1i * data2d(2);
        end
        [evals, order] = sort(evals, 'ComparisonMethod', 'abs');
		%fprintf('evals_first = %0.16e + i*%0.16e\n', real(evals(1)), imag(evals(1)));
        evals = evals(2:end);
        
        s_id = (seed_id - 1) * (N2 - 1) + 1;
        f_id = seed_id * (N2 - 1);
        all_evals(s_id : f_id) = evals;      
    end
    
    suffix = sprintf('ns(%d)_numSeeds(%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)', ...
        ns, ...
        num_seeds, ...
        diss_type, ...
        diss_phase, ...
        diss_gamma, ...
        W, ...
        U, ...
        J);
    
    pdf2d.x_bin_s = min(real(all_evals));
    pdf2d.x_bin_f = max(real(all_evals));
    pdf2d.y_bin_s = min(imag(all_evals));
    pdf2d.y_bin_f = max(imag(all_evals));
    pdf2d = oqs_pdf_2d_setup(pdf2d);
    data2d = horzcat(real(all_evals), imag(all_evals));
    pdf2d = oqs_pdf_2d_update(pdf2d, data2d);
    pdf2d = oqs_pdf_2d_release(pdf2d);
    fig = oqs_pdf_2d_plot(pdf2d);
    fn_fig = sprintf('%s/lindbladian_evals_%s', figures_path, suffix);
    oqs_save_fig(fig, fn_fig) 
end

