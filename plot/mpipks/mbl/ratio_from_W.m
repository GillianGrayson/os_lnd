clear all;
addpath('../../../source/matlab/lib')

path = '/data/condmat/ivanchen/yusipov/os_lnd/mbl';
figures_path = '/home/ivanchen/yusipov/os_lnd/figures/mbl';

num_T = 500;
step = 0.01;

N = 12;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

Ws = linspace(0.2, 20, 100)';
U = 1.0;
J = 1.0;

seeds = linspace(1, 100, 100)';

num_Ws = size(Ws, 1);
num_seeds = size(seeds, 1);
ratios = zeros(num_Ws, num_seeds);

ratio_pdf.xs = Ws;
ratio_pdf.x_num_points = num_Ws;
ratio_pdf.y_num_bins = num_Ws;
ratio_pdf.x_label = '$W$';
ratio_pdf.y_label = '$r$';

ratio_std_shade.xs = Ws;
ratio_std_shade.alpha = 0.35;
ratio_std_shade.legend = sprintf('N=%d', N);
ratio_std_shade.x_label = '$W$';
ratio_std_shade.y_label = '$r$';

for W_id = 1:num_Ws
    
    W = Ws(W_id);
    fprintf('W = %0.16e\n', W);
    
    for seed_id = 1:num_seeds
	
        seed = seeds(seed_id);
		
        suffix = sprintf('times(%0.2e_%0.2e)_ns(%d)_seed(%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)', ...
            0, ...
			num_T, ...
			N, ...
            seed, ...
            diss_type, ...
            diss_phase, ...
            diss_gamma, ...
            W, ...
            U, ...
            J);
        
        fn = sprintf('%s/odeint_%d_%0.2e/ns_%d/diss_%d_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/seed_%d/ratios_%s.txt', ...
            path, ...
            num_T, ...
            step, ...
            N, ...
            diss_type, ...
            diss_phase, ...
            diss_gamma, ...
            W, ...
            U, ...
            J, ...
            seed, ...
            suffix);
        
        data = importdata(fn);
        ratios(W_id, seed_id) = data(end);       
    end
end

suffix = sprintf('ns(%d)_numSeeds(%d)_diss(%d_%0.4f_%0.4f)_prm(var_%0.4f_%0.4f)', ...
    N, ...
    num_seeds, ...
    diss_type, ...
    diss_phase, ...
    diss_gamma, ...
    U, ...
    J);

ratio_pdf.y_bin_s = min(ratios, [], 'all');
ratio_pdf.y_bin_f = max(ratios, [], 'all');
ratio_pdf = oqs_pdf_2d_lead_by_x_setup(ratio_pdf);
ratio_pdf = oqs_pdf_2d_lead_by_x_update(ratio_pdf, ratios);
ratio_pdf = oqs_pdf_2d_lead_by_x_release(ratio_pdf);
fig = oqs_pdf_2d_lead_by_x_plot(ratio_pdf);
fn_fig = sprintf('%s/ratio_pdf_%s', figures_path, suffix);
oqs_save_fig(fig, fn_fig)

ratio_std_shade.data = ratios;
fig = oqs_plot_std_shade(ratio_std_shade);
fn_fig = sprintf('%s/ratio_std_shade_%s', figures_path, suffix);
oqs_save_fig(fig, fn_fig)
