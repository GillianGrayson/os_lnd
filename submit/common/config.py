def get_regular_global_config(system, task, save_rho='true', name_precision=4):
    config_list = []
    config_list.append('[global]')
    config_list.append('run_type = regular')
    config_list.append('logger_type = console')
    config_list.append('silent = false')
    config_list.append('system = ' + system)
    config_list.append('task = ' + task)
    config_list.append('debug_dump = false')
    config_list.append('save_hamiltonians = false')
    config_list.append('save_dissipators = false')
    config_list.append('save_f_basis = false')
    config_list.append('save_lindbladians = false')
    config_list.append('save_lindbladian_evals = true')
    config_list.append('save_rho = ' + save_rho)
    config_list.append('save_rho_evals = true')
    config_list.append('save_run_times = true')
    config_list.append('save_non_zeros_part = true')
    config_list.append('save_mem_info = true')
    config_list.append('name_precision = ' + str(name_precision))
    config_list.append('save_precision = 16')
    return config_list

def get_serial_global_config(system, task, start, shift, num, name_precision=4):
    config_list = []
    config_list.append('[global]')
    config_list.append('run_type = serial')
    config_list.append('serial_start = ' + str(start))
    config_list.append('serial_shift = ' + str(shift))
    config_list.append('serial_num = ' + str(num))
    config_list.append('logger_type = console')
    config_list.append('silent = false')
    config_list.append('system = ' + system)
    config_list.append('task = ' + task)
    config_list.append('debug_dump = false')
    config_list.append('save_hamiltonians = false')
    config_list.append('save_dissipators = false')
    config_list.append('save_f_basis = false')
    config_list.append('save_lindbladians = false')
    config_list.append('save_lindbladian_evals = false')
    config_list.append('save_rho = false')
    config_list.append('save_rho_evals = false')
    config_list.append('save_run_times = false')
    config_list.append('save_non_zeros_part = false')
    config_list.append('save_mem_info = false')
    config_list.append('name_precision = ' + str(name_precision))
    config_list.append('save_precision = 16')
    return config_list


def get_odeint_config(
        step,
        num_obser_periods,
        num_trans_periods,
        current_num_obser_periods,
        current_num_obser_time_points,
        is_continue,
        continue_path,
        dump_progress,
        dump_last_time
):
    config_list = []
    config_list.append('[odeint]')
    config_list.append('start_state_type = 1')
    config_list.append('start_state_id = 0')
    config_list.append('step = ' + str(step))
    config_list.append('dump_type = linear')
    config_list.append('num_trans_periods = ' + str(num_trans_periods))
    config_list.append('num_obser_periods = ' + str(num_obser_periods))
    config_list.append('current_num_obser_periods = ' + str(current_num_obser_periods))
    config_list.append('current_num_obser_time_points = ' + str(current_num_obser_time_points))
    config_list.append('dump_progress = ' + str(dump_progress))
    config_list.append('dump_last_time = ' + str(dump_last_time))
    config_list.append('continue = ' + str(is_continue))
    config_list.append('continue_path = ' + str(continue_path))
    return config_list


def get_smallest_eigen_vector_config(max_num_iterations, tolerance):
    config_list = []
    config_list.append('[smallest_eigen_vector]')
    config_list.append('max_num_iterations = ' + str(max_num_iterations))
    config_list.append('tolerance = ' + str(tolerance))
    return config_list


def get_all_evals_config(max_num_iterations, tolerance):
    config_list = []
    config_list.append('[all_evals]')
    config_list.append('max_num_iterations = ' + str(max_num_iterations))
    config_list.append('tolerance = ' + str(tolerance))
    return config_list