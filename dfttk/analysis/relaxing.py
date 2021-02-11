import numpy as np

def transmat(c1, c2):

    """ Get the transformation matrix between c1 and c2.

    Defined such that:
    Lattice(transmat(initial.lattice.matrix, final.lattice.matrix)@initial.lattice.matrix).scale(final.volume).matrix == final.lattice.matrix
    """
    c1 = c1/(np.linalg.det(c1)**(1.0/3.0))
    c2 = c2/(np.linalg.det(c2)**(1.0/3.0))
    return np.linalg.inv(c1)@c2

def get_non_isotropic_strain(c1, c2):
    """Based on ATAT's implementation for checkrelax, strain required to transform cell c1 to c2, neglecting isotropic scaling or rotations.
    
    c1 and c2 should both be structure.lattice.matrix instances
    """
    return np.linalg.norm(transmat(c1, c2) - np.eye(3))


def get_bond_distance_change(initial_structure, final_structure):
    """Square root of sum of square distance matrices (controlling for cell shape change)"""
    # Divides by 2 because the bonds are double counted in the symmetric distance matrix
    bond = np.linalg.norm((initial_structure.distance_matrix - final_structure.distance_matrix)/2)/len(initial_structure)
    return bond
