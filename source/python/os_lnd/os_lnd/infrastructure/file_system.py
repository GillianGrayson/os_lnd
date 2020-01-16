import socket
import getpass

def get_input_path(params_path):

    host_name = socket.gethostname()
    path = ''
    if host_name == 'MSI':
        path = 'D:/Work/os_lnd/source/cpp/os_lnd/os_lnd'
    elif host_name == 'DESKTOP-K9VO2TI':
        path = 'E:/Work/os_lnd/source/cpp/os_lnd/os_lnd'
    elif host_name == 'master' or host_name[0:4] == 'node':
        user = getpass.getuser()
        if user == 'ivanchen':
            path = '/data3/ivanchen/yusipov/os_lnd/' + params_path

    return path
