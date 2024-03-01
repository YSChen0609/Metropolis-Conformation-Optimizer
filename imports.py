# imports.py

## Just some tools. 

import time
import numpy as np
import itertools

def log_method(func):
    def wrapper(*args, **kwargs):
        print(f"Starting {func.__name__}...")
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"Finished {func.__name__}. Time taken: {end_time - start_time} seconds.")
        return result
    return wrapper


def angles_between(vec1, vec2):
    """
    Get the angles (degrees) between vec1 and vec2.
    Note that the vectors are ndarrays.
    """
    return np.degrees(
                np.arccos(
                    np.dot(vec1/np.linalg.norm(vec1), vec2/np.linalg.norm(vec2))
                )
            )

def get_pair_idx(conformation_info=np.zeros((5,3))):
    """
    Get the combination (n choose 2) of indices.
    Note that we skip the first row (which serves as the origin atom).

    - Assume the first row being the atom on the origin.
    - conformation_info: ndarray with [atom_name(str): skipped], x, y, z.
    - default shape = (5,3)
    """
    idx = range(1,conformation_info.shape[0]) # skip the first row (origin atom)

    return tuple(itertools.combinations(idx,2))
