import configparser
from os_lnd.model.model import Model
from os_lnd.model.context import Context
from os_lnd.infrastructure.file_system import get_input_path

params_path = ''
ini_fn = get_input_path(params_path) + '/config.ini'
ini = configparser.ConfigParser()
ini.read(ini_fn)
model = Model(ini)
context = Context(model)
context.run(model)