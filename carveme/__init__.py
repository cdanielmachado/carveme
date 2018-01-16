import os
from ConfigParser import ConfigParser
from framed import set_default_solver, set_default_parameter, Parameter

__version__ = '1.1.0'

project_dir = os.path.abspath(os.path.dirname(__file__)) + '/'

config = ConfigParser()
config.read(project_dir + 'config.cfg')

set_default_solver(config.get('solver', 'default_solver'))
set_default_parameter(Parameter.FEASIBILITY_TOL, config.getfloat('solver', 'feas_tol'))
set_default_parameter(Parameter.OPTIMALITY_TOL, config.getfloat('solver', 'opt_tol'))
set_default_parameter(Parameter.INT_FEASIBILITY_TOL, config.getfloat('solver', 'int_feas_tol'))