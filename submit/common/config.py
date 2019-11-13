def get_global_config(system, task):
    config_list = []
    config_list.append('[global]')
    config_list.append('logger_type = console')
    config_list.append('silent = false')
    config_list.append('system = ' + system)
    config_list.append('task = ' + task)
    config_list.append('debug_dump = true')
    config_list.append('name_precision = 4')
    config_list.append('save_precision = 16')
    return config_list

def get_odeint_config(step, start_observed_period, finish_observed_period):

    num_time_points = finish_observed_period - start_observed_period + 1

    config_list = []
    config_list.append('[odeint]')
    config_list.append('start_state_id = 0')
    config_list.append('step = ' + str(step))
    config_list.append('dump_type = linear')
    config_list.append('start_observed_period = ' + str(start_observed_period))
    config_list.append('finish_observed_period = ' + str(finish_observed_period))
    config_list.append('num_time_points = ' + str(num_time_points))
    config_list.append('dump_progress = true')
    return config_list
