# config.py

## Defines the parameters used to initialize MCMC objects (and for experiments).

from collections import namedtuple

Config = namedtuple('Config', 
                    ['pdb_file_path',
                     'Kr',                      # spring_constant for length
                     'r0',                      # equilibrium length
                     'K_theta',                 # spring_constant for angle
                     'theta0',                  # equilibrium angle
                     'kT_list',                 # Initial kT
                     'perturb_scale_gaussian',  # The std for gaussian noise
                     'iter_rounds',             # iteration rounds for metrpolis algo
                     'record_interval'          # number of record every n iterations
                    ])