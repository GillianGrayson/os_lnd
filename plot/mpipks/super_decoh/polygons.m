clear all;
addpath('../../../source/matlab/lib')

%ps = vertcat([0.01]', [0.1:0.02:0.3]', [0.5, 0.8, 1.0]');
ps = [0.1];
type = 'source';

reshufle_type = 1;
G_type = 0;
N = 150;
seeds = linspace(1, 100, 100)';

num_hits = N/3;

if type == 'source'
	theory_data = importdata(sprintf('borderline_classical.dat'));
else
	theory_data = importdata(sprintf('borderline_%s.dat', type'));
end

xor_areas = zeros(size(ps, 1), 1);
xor_areas_normed = zeros(size(ps, 1), 1);

for p_id = 1:size(ps, 1)

	p = ps(p_id);

	fprintf('p = %0.16e\n', p);
    
    pdf2d.x_num_bins = 201;
    pdf2d.y_num_bins = 201;
    pdf2d.x_label = '$Re(\lambda)$';
    pdf2d.y_label = '$Im(\lambda)$';
    
    path = '/data/condmat/ivanchen/yusipov/os_lnd/super_decoh/eigen_dense';
    figures_path = '/home/ivanchen/yusipov/os_lnd/figures/super_decoh';
    
    N2 = N * N;
    
    all_evals = zeros((N2 - 1) * size(seeds, 1), 1);
    
    for seed_id = 1:size(seeds, 1)
        
        seed = seeds(seed_id);
        fprintf('seed = %d\n', seed);
        
        suffix = sprintf('reshuffle(%d)_G(%d)_N(%d)_p(%0.10f)_seed(%d)', reshufle_type, G_type, N, p, seed);
        
        evals = zeros(N2, 1);
        fn_cpp = sprintf('%s/G_type_%d/reshuffle_type_%d/N_%d/p_%0.10f/seed_%d/lindbladian_evals_%s.txt', path, G_type, reshufle_type, N, p, seed, suffix);
        if ~isfile(fn_cpp)
			fprintf('Warning! seed = %d\n', seed);
		end
		cpp_data = importdata(fn_cpp);
        for str_id = 1:N2
            str = string(cpp_data(str_id));
            data2d = sscanf(str, '(%e,%e)', 2);
            evals(str_id) = data2d(1) + 1i * data2d(2);
        end
        
        [evals, order] = sort(evals, 'ComparisonMethod', 'abs');
        evals = evals(2:end);
		if strcmp(type, 'classical')
			evals = N * sqrt(N) * (real(evals) + 1) + 1i * N * sqrt(N) * imag(evals);
        elseif strcmp(type, 'quantum')
			evals = N / p * (real(evals) + 1) + 1i * N / p * imag(evals);
        elseif strcmp(type, 'source')
			evals = evals;
		else
            error('Wrong type');
        end
        
        s_id = (seed_id - 1) * (N2 - 1) + 1;
        f_id = seed_id * (N2 - 1);
        all_evals(s_id : f_id) = evals;
    end
     
    suffix = sprintf('numHits(%d)_type(%s)_reshuffle(%d)_N(%d)_p(%0.10f)_numSeeds(%d)', num_hits, type, reshufle_type, N, p, size(seeds, 1));
    
    pdf2d.x_bin_s = min(real(all_evals)) - 1e-16;
    pdf2d.x_bin_f = max(real(all_evals)) + 1e-16;
    pdf2d.y_bin_s = min(imag(all_evals)) - 1e-16;
    pdf2d.y_bin_f = max(imag(all_evals)) + 1e-16;
    pdf2d = oqs_pdf_2d_setup(pdf2d);
    data2d = horzcat(real(all_evals), imag(all_evals));
    pdf2d = oqs_pdf_2d_update(pdf2d, data2d);
    pdf2d = oqs_pdf_2d_release(pdf2d);
    fig = oqs_pdf_2d_plot(pdf2d);
    fn_fig = sprintf('%s/lindbladian_evals_%s', figures_path, suffix);
    
    xt1 = theory_data(:, 1);
    yt1 = theory_data(:, 2);
    xt2 = flip(xt1);
    yt2 = -flip(yt1);
    xt = vertcat(xt1, xt2);
    yt = vertcat(yt1, yt2);
    
    pgont = polyshape(xt, yt);
    
    hold all;
    plot(pgont);
	
	limit_pdf = num_hits / (size(all_evals, 1) * pdf2d.x_bin_shift * pdf2d.y_bin_shift)
	
    xr1 = pdf2d.x_bin_centers';
    yr1 = zeros(size(pdf2d.x_bin_centers, 2), 1);
    for x_id = 1:size(pdf2d.x_bin_centers, 2)
        is_found = 0;
        for y_id = 1:size(pdf2d.y_bin_centers, 2)
            if pdf2d.pdf(x_id, y_id) > limit_pdf
                yr1(x_id) = pdf2d.y_bin_centers(y_id);
                is_found = 1;
                break;
            end
            if is_found == 0
                yr1(x_id) = 0;
            end
        end
    end
    
    yr1 = smooth(yr1 * -1.0);
    xr2 = flip(xr1);
    yr2 = -flip(yr1);
    
    xr = vertcat(xr1, xr2);
    yr = vertcat(yr1, yr2);
    pgonr = polyshape(xr, yr);
	
	contour = horzcat(xr, yr);
    
	fn_txt = sprintf('%s/contour_%s.txt', figures_path, suffix);
	%save(fn_txt, 'contour', '-double', '-tab');
    %dlmwrite(sprintf('%s/contour_%s.txt', figures_path, suffix), contour);
	writematrix(contour, fn_txt, 'Delimiter','tab')
	
    hold all;
    plot(pgonr);
    
    polyout = xor(pgont,pgonr);
    
    hold all;
    plot(polyout);
	
	polyunion = union(pgont,pgonr);

    xor_areas(p_id) = area(polyout);
	xor_areas_normed(p_id) = xor_areas(p_id) / area(polyunion);
	
	fprintf('xor_areas_normed = %0.16e\n', xor_areas_normed(p_id));
    
    oqs_save_fig(fig, fn_fig);
    
end

suffix = sprintf('numHits(%d)_type(%s)_reshuffle(%d)_N(%d)_p(var)_numSeeds(%d)', num_hits, type, reshufle_type, N, size(seeds, 1));

fig = figure;
plot(ps, xor_areas, 'LineWidth', 2)
set(gca, 'FontSize', 30);
xlabel('p', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('D', 'Interpreter', 'latex');
set(gca,'XScale','log');
set(gca,'YScale','log');
hold all;
fn_fig = sprintf('%s/diff_from_p_%s', figures_path, suffix);
oqs_save_fig(fig, fn_fig)

fig = figure;
plot(ps, xor_areas_normed, 'LineWidth', 2)
set(gca, 'FontSize', 30);
xlabel('p', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('D', 'Interpreter', 'latex');
set(gca,'XScale','log');
set(gca,'YScale','log');
hold all;
fn_fig = sprintf('%s/diff_normed_from_p_%s', figures_path, suffix);
oqs_save_fig(fig, fn_fig)
