# MCMC.py

import pandas as pd
import numpy as np

from imports import *

class MCMC():
    """
    Implement the Metropolis algorithm to optimize (find) the conformation with a starting conformation.
    We only perturbate the vectors (x,y,z) of the existing atoms, no atom insertion/deletion is involved.
    """
    def __init__(self, config):
        self.config = config
        self.starting_conformation = self.readPDB(self.config.pdb_file_path)
        self.angle_pair_idx =  get_pair_idx(self.starting_conformation)

    def metropolis_algo(self, kT):
        """
        Implement Metropolis Algorithm.
        """
        # Initialization
        E_i    = self.getEnergy(self.starting_conformation)
        E_min  = self.getEnergy(self.starting_conformation)
        X_i    = self.starting_conformation.copy() # this sample
        X_min  = self.starting_conformation.copy()
        result = [{
                    'round': 1,
                    'E_min': E_min, # lowest energy so far
                    'kT'   : kT
                  }]

        # Execute the Metropolis Algo. for certain iterations
        for i in range(1, self.config.iter_rounds+1):
            E_i, X_i, Accept = self.metropolis_step(E_i, X_i, kT)

            if Accept:
                E_min = E_i
                X_min = X_i

            if i % self.config.record_interval == 0: #record every 'record_interval' iterations
                result.append(
                    {
                        'round': i,
                        'E_min': E_min, # lowest energy so far
                        'kT'   : kT
                    }
                )

        return E_min, X_min, pd.DataFrame.from_dict(result), kT


    def metropolis_step(self, E_i, X_i, kT=0.6):
        """
        The subroutine for Metropolis Algorithm (self.metropolis_algo):
            1. Generate a perturbation (trial sample)
            2. Accept with prob = min(1, exp(-1/kT * (E_curr - E_i)))
        """
        Accept = True # return an indicator of whether accept or not

        # Generate a perturbation (trial sample)
        trial_sample = self.perturbate(X_i)
        E_trial = self.getEnergy(trial_sample)

        # Energy is lower, accept
        if  E_trial < E_i:
            # accept and return the current info
            return E_trial, trial_sample, Accept

        # Energy is Not lower, accept with a probability
        else:
            acc_prob = np.exp(-1/kT * (E_trial - E_i))

            if np.random.uniform(0, 1) <= acc_prob:
                # accept
                return E_trial, trial_sample, Accept
            else:
                # reject
                return E_i, X_i, ~Accept


    def perturbate(self, molecular_info):
        """
        Make a perturbation and generate a trial sample.
        Here we are using a simple gaussian noise.
        """
        noise = np.concatenate((np.zeros((1,molecular_info.shape[1]), dtype=float),np.random.normal(loc=0.0, scale=self.config.perturb_scale_gaussian, size=molecular_info[1:].shape)))
        trial_sample = molecular_info + noise

        return trial_sample

    def getEnergy(self, molecular_info):
        """
        Calculate the (simplified) molecular mechanics force field given the molecular conformation table.
            - E = sum(all bonds) 1/2 Kr * (r - r0)^2 + sum(all angles) 1/2 Kθ * (θ - θ0)^2

        Note that here we fix a specific atom to the origin, since we are using relative molecular coordinates.
        """
        E_bond  = (1/2) * self.config.Kr *\
                    np.apply_along_axis(lambda row: (np.linalg.norm(row) - self.config.r0)**2,
                                        axis=1, arr=molecular_info).sum()

        E_angle = (1/2) * self.config.K_theta *\
                    np.array([(angles_between(*molecular_info[pair,:]) - self.config.theta0)**2
                              for pair in self.angle_pair_idx]).sum()

        return E_bond + E_angle


    def readPDB(self, pdb_file_path):
        """
        Read and parse the starting conformation file (.pdb).
        """
        segData = []
        with open(pdb_file_path, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith('ATOM'):
                    segData.append({
#                         'atom_name' : line[12:16].strip(),
                        'x_coord'   : float(line[30:38]),
                        'y_coord'   : float(line[38:46]),
                        'z_coord'   : float(line[46:54])
                    })

            return pd.DataFrame.from_dict(segData).to_numpy(dtype=float)


