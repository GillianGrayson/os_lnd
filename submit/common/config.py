def get_global_config(system, task, name_precision=4):
    config_list = []
    config_list.append('[global]')
    config_list.append('logger_type = console')
    config_list.append('silent = false')
    config_list.append('system = ' + system)
    config_list.append('task = ' + task)
    config_list.append('debug_dump = false')
    config_list.append('save_hamiltonians = false')
    config_list.append('save_dissipators = false')
    config_list.append('save_lindbladians = false')
    config_list.append('save_f_basis = false')
    config_list.append('name_precision = ' + str(name_precision))
    config_list.append('save_precision = 16')
    return config_list

def get_odeint_config(step, total_num_periods, current_num_periods, current_num_time_points, is_continue, continue_path):

    config_list = []
    config_list.append('[odeint]')
    config_list.append('start_state_id = 0')
    config_list.append('step = ' + str(step))
    config_list.append('dump_type = linear')
    config_list.append('total_num_periods = ' + str(total_num_periods))
    config_list.append('current_num_periods = ' + str(current_num_periods))
    config_list.append('current_num_time_points = ' + str(current_num_time_points))
    config_list.append('dump_progress = true')
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