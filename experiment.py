# experiment.py

## Summarize and Plot the experiment result given a experiment configuration.

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd


def experiment(config):
    def plot_trajectory(records):
        """
        Plot the trajectory, showing the lowest energy so far for every certain iterations.
        """
        sns.lineplot(data=records, x='round', y='E_min', hue='kT')

        plt.title('Lowest energy so far for every 1000 iterations')
        plt.ylabel('Energy (kcal)')

        plt.legend(title='kT (kcal/mol)')
        plt.show()

    def get_summary_table(experiment_result):
        """
        Return a summary table of the experiment result.
        Format: kT, E_min, [bond_lengths], [angles] #TODO: angles not general?

        Note that this is specify for methane experiment, need to customize to your task.
        """
        summary_rows = []
        angle_pair_idx = get_pair_idx()
        
        for result in experiment_result:
            bond_lengths = np.apply_along_axis(lambda row: np.linalg.norm(row),
                                 axis=1, arr=result[1][1:,:])
            angles    = np.array([angles_between(*result[1][pair,:])
                            for pair in angle_pair_idx])
        
            summary_row = {
                'kT': result[3],
                'E_min': result[0],
                **{f'bond_{i}': bond_length for i, bond_length in enumerate(bond_lengths, start=1)},
                **{f'angle_{i}': angle for i, angle in enumerate(angles, start=1)}
                }
        
            summary_rows.append(summary_row)

        return pd.DataFrame(summary_rows)

    @log_method
    def do_experiment(kT_list):
    """
    Experiment different kTs using Metropolis Algorithm.

    - metropolis_algo return format: E_min, X_min, records(DataFrame), kT
    - records format: round, E_min, kT
    """
        experiment_result = []

        for kT in kT_list:
            m = MCMC(config)
            experiment_result.append((m.metropolis_algo(kT)))

        global_best_result = min(experiment_result, key=lambda x: x[0])

        records = pd.concat([x[2] for x in experiment_result])

        summary_table = get_summary_table(experiment_result)

        # display the experiment results
        print(f"Lowest Energy (among all kTs): {global_best_result[0]}")
        display(summary_table)
        plot_trajectory(records)

    # Experiment starts here
    do_experiment(config.kT_list)