clear all;

N = 5;

imag_lim = 5e-3;

num_combos = 4;
num_states = 2^N;
M = num_states.^2 - 1;
all_evals = zeros(M * num_combos, 1);
zs_all = zeros(M * num_combos, 1);
s_id = 0;
f_id = 0;
for combo = 1:num_combos
    
    fprintf('seed = %d\n', combo);
    if (combo == 1)
        tau = 1;
        k = 1;
    elseif (combo == 2)
        tau = 1;
        k = -1;
    elseif (combo == 3)
        tau = -1;
        k = 1;
    else
        tau = -1;
        k = -1;
    end
    
    con = zeros(4);
    con(1, 1) = tau;
    con(4, 4) = k;
    con(2, 3) = 1;
    con(3, 2) = 1;
    
    dissipators = cell(N, 1);
    
    
    % =====================================================================
    %                       Left dissipator (1,2)
    % =====================================================================
    left_diss = con;
    for spin_id = 3:N
        left_diss = kron(left_diss, eye(2));
    end
    dissipators{1} = left_diss;
    
    % =====================================================================
    %                      Right dissipator (N-1, N)
    % =====================================================================
    right_diss = eye(2);
    for spin_id = 2:N-2
        right_diss = kron(right_diss, eye(2));
    end
    right_diss = kron(right_diss, con);
    dissipators{N-1} = right_diss;
    
    % =====================================================================
    %                  All middle dissipators (N-1, N)
    % =====================================================================
    for diss_id = 2:N-2
        
        curr_diss = eye(2);
        
        left_s_id = 2;
        left_f_id = diss_id - 1;
        for spin_id = left_s_id : left_f_id
            curr_diss = kron(curr_diss, eye(2));
        end
        
        curr_diss = kron(curr_diss, con);
        
        right_s_id = diss_id + 2;
        right_f_id = N;
        for spin_id = right_s_id : right_f_id
            curr_diss = kron(curr_diss, eye(2));
        end
        
        dissipators{diss_id} = curr_diss;
    end
    
    % =====================================================================
    %                       Border dissipator (N, 1)
    % =====================================================================
    border_diss = zeros(num_states);
    for state_id = 1:num_states
        state_int = state_id - 1;
        state_bits = de2bi(state_int, N);
        state_bits = flip(state_bits);
        
        if (state_bits(1) == state_bits(end))
            if state_bits(1) == 0
                border_diss(state_id, state_id) = tau;
            else
                border_diss(state_id, state_id) = k;
            end
        else
            inv_bits = state_bits;
            inv_bits(1) = state_bits(end);
            inv_bits(end) = state_bits(1);
            inv_bits = flip(inv_bits);
            row_id = bi2de(inv_bits) + 1;
            border_diss(row_id, state_id) = 1;
        end
    end
    dissipators{N} = border_diss;
    
    % =====================================================================
    %      Check the first dissipator which are constructed this way
    % =====================================================================
    first_diss = zeros(num_states);
    for state_id = 1:num_states
        state_int = state_id - 1;
        state_bits = de2bi(state_int, N);
        state_bits = flip(state_bits);
        
        if (state_bits(1) == state_bits(2))
            if state_bits(1) == 0
                first_diss(state_id, state_id) = tau;
            else
                first_diss(state_id, state_id) = k;
            end
        else
            inv_bits = state_bits;
            inv_bits(1) = state_bits(2);
            inv_bits(2) = state_bits(1);
            inv_bits = flip(inv_bits);
            row_id = bi2de(inv_bits) + 1;
            first_diss(row_id, state_id) = 1;
        end
    end
    fprintf('first dissipator difference = %0.4e\n', norm(dissipators{1} - first_diss));
    
    % =====================================================================
    %                            Lindbladian
    % =====================================================================
    Lind = zeros(num_states^2);
    for diss_id = 1:size(dissipators, 1)
        tmp = 0.5 * (2 * kron(eye(num_states), dissipators{diss_id}) * ...
            kron(transpose(dissipators{diss_id}'), eye(num_states)) - ...
            kron(transpose(dissipators{diss_id}' * dissipators{diss_id}), eye(num_states)) - ...
            kron(eye(num_states), dissipators{diss_id}' * dissipators{diss_id}));
        Lind = Lind + tmp;
    end
    
    
    % =====================================================================
    %                            Evals routines
    % =====================================================================
    evals = eig(Lind);
    evals = sort(evals,'ComparisonMethod','abs');
    evals = evals(2:end);
    all_evals((combo - 1) * M + 1 : combo * M) = evals;
    
    evals_filtered = evals(abs(imag(evals)) >= imag_lim);
    num_passed_evals = size(evals_filtered, 1);
    fprintf('num_passed_evals = %d\n', num_passed_evals);
    
    s_id = f_id + 1;
    f_id = s_id + num_passed_evals - 1;
    
    zs = zeros(num_passed_evals, 1);
    neighbours_2D = horzcat(real(evals_filtered), imag(evals_filtered));
    for z_id = 1 : num_passed_evals
        target_2D = horzcat(real(evals_filtered(z_id)), imag(evals_filtered(z_id)));
        target = evals_filtered(z_id);
        distances = sqrt((neighbours_2D(:, 1) - target_2D(:, 1)).^2 + (neighbours_2D(:, 2) - target_2D(:, 2)).^2);
        [min_distances, order] = mink(distances, 3);
        nn = evals_filtered(order(2));
        nnn = evals_filtered(order(3));
        zs(z_id) = (nn - target) / (nnn - target);
    end
    
    num_zero_zs = sum(abs(zs) < 1e-4);
    fprintf('num_zero_zs = %d\n', num_zero_zs);
    
    zs_all(s_id : f_id) = zs;
end

zs_all = zs_all(1 : f_id);

total_num_passed_evals = f_id

pdf2dzs.x_num_bins = 101;
pdf2dzs.y_num_bins = 101;
pdf2dzs.x_label = '$Re(z)$';
pdf2dzs.y_label = '$Im(z)$';

pdf2dzs.x_bin_s = -1;
pdf2dzs.x_bin_f = 1;
pdf2dzs.y_bin_s = -1;
pdf2dzs.y_bin_f = 1;
pdf2dzs = oqs_pdf_2d_setup(pdf2dzs);
data2d = horzcat(real(zs_all), imag(zs_all));
pdf2dzs = oqs_pdf_2d_update(pdf2dzs, data2d);
pdf2dzs = oqs_pdf_2d_release(pdf2dzs);
fig = oqs_pdf_2d_plot(pdf2dzs);
fn_fig = sprintf('zs_N(%d)_num_seeds(%d)', N, num_combos);
oqs_save_fig(fig, fn_fig)

pdf2d.x_num_bins = 201;
pdf2d.y_num_bins = 201;
pdf2d.x_label = '$Re(\lambda)$';
pdf2d.y_label = '$Im(\lambda)$';

pdf2d.x_bin_s = min(real(all_evals));
pdf2d.x_bin_f = max(real(all_evals));
pdf2d.y_bin_s = min(imag(all_evals));
pdf2d.y_bin_f = max(imag(all_evals));

pdf2d = oqs_pdf_2d_setup(pdf2d);
data2d = horzcat(real(all_evals), imag(all_evals));
pdf2d = oqs_pdf_2d_update(pdf2d, data2d);
pdf2d = oqs_pdf_2d_release(pdf2d);
fig = oqs_pdf_2d_plot(pdf2d);
fn_fig = sprintf('Levals_N(%d)_num_seeds(%d)', N, num_combos);
oqs_save_fig(fig, fn_fig)

