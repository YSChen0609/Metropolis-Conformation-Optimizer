# __main__.py

from config import config
from experiment import experiment


def main():
    # Config for experimenting optimization for methane
    
    # kT_list is decided based on the following:
    # The melting point of methane is approximately 90 Kelvin. (kT=0.179)
    # The boiling point of methane is approximately 112 Kelvin (kT=0.2224).
    # The auto-ignition point of methane is approximately 813 Kelvin (kT=1.615).
    
    config = Config(pdb_file_path='./data/methane_start.pdb',
                    Kr=367, r0=1.08, K_theta=35, theta0=109.5, 
                    kT_list=[0.179,0.2224,0.6,1.615],
                    perturb_scale_gaussian=.75,
                    iter_rounds = 100000, record_interval=1000)

    experiment(config)

if __name__ == "__main__":
    main()