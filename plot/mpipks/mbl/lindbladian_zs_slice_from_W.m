clear all;
addpath('../../../source/matlab/lib')

task = 'eigen_dense';

path = sprintf('/data/condmat/ivanchen/yusipov/os_lnd/mbl/%s', task);

figures_path = '/home/ivanchen/yusipov/os_lnd/figures/mbl';

without_line = 0;

ns = 8;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

Ws = linspace(0.2, 20, 100)';
U = 1.0;
J = 1.0;

global_seeds = 250;
target_seeds = 250;

seeds = linspace(1, target_seeds, target_seeds)';

num_Ws = size(Ws, 1);
num_seeds = size(seeds, 1);

zs_slice_pdf.xs = Ws;
zs_slice_pdf.x_num_points = num_Ws;
zs_slice_pdf.y_num_bins = num_Ws;
zs_slice_pdf.x_label = '$W$';
zs_slice_pdf.y_label = '$z slice$';
zs_slice_pdf.y_bin_s = -1;
zs_slice_pdf.y_bin_f = +1;
zs_slice_pdf = oqs_pdf_2d_lead_by_x_setup(zs_slice_pdf);

N = (nchoosek(ns, ns/2));
N2 = N^2;

for W_id = 1:num_Ws
    
    W = Ws(W_id);
    fprintf('W = %0.16e\n', W);

    zs_slice_all = zeros(N2 * size(seeds, 1), 1);
       
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
    
    suffix = sprintf('ns(%d)_seeds(%d_%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)', ...
        ns, ...
        1, ...
        1, ...
        global_seeds, ...
        diss_type, ...
        diss_phase, ...
        diss_gamma, ...
        W, ...
        U, ...
        J);
    
    fn = sprintf('%s/lindbladian_evals_%s.txt', aggr_path, suffix);

    lind_evals_all = importdata(fn);
    
    curr_fin_id = 0;
    curr_fin_slice_id = 0;
    
    for seed_id = 1:num_seeds
	
        seed = seeds(seed_id);
		fprintf('seed = %d\n', seed);
		
        evals = lind_evals_all((seed_id - 1) * N2 + 1 : seed_id * N2, 1) + 1i * lind_evals_all((seed_id - 1) * N2 + 1 : seed_id * N2, 2);
        [evals, order] = sort(evals, 'ComparisonMethod', 'abs');
        evals = evals(2:end);
        
        if without_line == 1
            evals = evals(abs(imag(evals(:))) > 1e-2); 
        end

        zs_slice = [];
        
        neighbours_2D = horzcat(real(evals), imag(evals));
        for z_id = 1 : size(evals, 1)
            target_2D = horzcat(real(evals(z_id)), imag(evals(z_id)));
            target = evals(z_id);
            distances = sqrt(sum(bsxfun(@minus, neighbours_2D, target_2D).^2, 2));
            [sorted_distances, order] = sort(distances);
            nn = evals(order(2));
            nnn = evals(order(3));
            z = (nn - target) / (nnn - target);
            
            if abs(imag(z)) < 5e-2
                zs_slice = vertcat(zs_slice, z);
            end
        end
        
        zs_slice_all(curr_fin_slice_id + 1 : curr_fin_slice_id + size(zs_slice, 1)) = zs_slice;
        curr_fin_slice_id = curr_fin_slice_id + size(zs_slice, 1);
        
    end

    zs_slice_all = zs_slice_all(1 : curr_fin_slice_id);
    
    zs_slice_pdf = oqs_pdf_2d_lead_by_x_update_slice(zs_slice_pdf, real(zs_slice_all), W_id);
    
    suffix = sprintf('ns(%d)_numSeeds(%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)', ...
        ns, ...
        num_seeds, ...
        diss_type, ...
        diss_phase, ...
        diss_gamma, ...
        W, ...
        U, ...
        J);
end

zs_slice_pdf = oqs_pdf_2d_lead_by_x_release(zs_slice_pdf, 0, 0);
fig = oqs_pdf_2d_lead_by_x_plot(zs_slice_pdf);
fn_fig = sprintf('%s/zs_slice_pdf_%s', figures_path, suffix);
oqs_save_fig(fig, fn_fig)

