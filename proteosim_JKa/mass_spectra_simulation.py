import proteosim_JKa as ps

amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}

def calculate_mol_mass(peptide_seq, amino_acid_mass_dalton):
    """
    Calculates the molecular mass of a peptide sequence and returns
    a dictionary containing the peptide sequence and its mass.

    Parameters
    ----------
    peptide_seq : str
        Peptide sequence consisting of one-letter amino acid codes.

    amino_acid_mass_dalton : dict
        Dictionary mapping amino acid one-letter codes to their masses in Dalton (Da).

    Returns
    -------
    dict
        Dictionary with the peptide sequence as key and the molecular
        mass as value.
    """
    mol_mass_dict = {}
    mol_mass = 0

    for amino_acid in peptide_seq:
        mol_mass += amino_acid_mass_dalton[amino_acid]

    mol_mass_dict[peptide_seq] = mol_mass

    return mol_mass_dict

def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
    Takes a str or list of str containing a aa sequence and returnes a dictionary of mass in dalton
    
    Parameters
    ----------
    peptides: str, list[str]
        List of or singular peptide sequences in one letter code
    amino_acid_mass_dict: dict
        Dictionary mapping amino acid one-letter codes to their masses in Dalton (Da).
    Returns
    -------
    dict
        Dictionary with the peptide sequence as key and the molecular mass as value. 
    """

    mol_mass_dict = {}

    # Falls nur ein Peptid übergeben wurde
    if isinstance(peptides, str):
        peptides = [peptides]

    for peptide_seq in peptides:
        mol_mass = 0

        for amino_acid in peptide_seq:
            mol_mass += amino_acid_mass_dict[amino_acid]

        mol_mass_dict[peptide_seq] = mol_mass

    return mol_mass_dict

def calculate_mz_collection(peptide_mass_map, charge, proton_mass=1.007):
    """
    Calculates the mass-to-charge ratio (m/z) for a collection of peptides
    from their molecular masses and returns the results as a dictionary.

    Parameters
    ----------
    peptide_mass_map : dict
        Dictionary mapping peptide sequences to their molecular masses in Dalton (Da)

    charge : int
        Charge state (z) used for the m/z calculation

    proton_mass : float, optional
        Mass of a proton in Dalton (Da). Default is 1.007 Da

    Returns
    -------
    dict
        Dictionary mapping peptide sequences to their m/z values for the specified charge state
    """

    dict_mz_peptides = {}
    for sequence, mass in peptide_mass_map.items():
        
        m_over_z = m_over_z = (mass + charge * proton_mass) / charge

        dict_mz_peptides[sequence] = {
            "mz" : m_over_z,  # hier potentiell noch
            "mass" : mass 
        }

    return dict_mz_peptides

import numpy as np
import matplotlib.pyplot as plt

def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42):
    """
    Generates a simple mass spectrum from a collection of m/z values by assigning random intensity values and visualizing the resulting peaks as a bar plot

    Parameters
    ----------
    mz_values : list[float]
        List of mass-to-charge (m/z) values to be displayed on the x-axis

    random_count_range : tuple[int, int], optional
        Minimum and maximum intensity values used to generate random peak heights; Default is (0, 30000)

    seed : int, optional
        Seed for the random number generator to ensure reproducible intensity values; Default is 42

    Returns
    -------
    plt.bar
        barplot of the generated spectrum
    """
    np.random.seed(seed)

    mz_values_list = [v["mz"] for v in mz_values.values()]

    intensities = [
        np.random.randint(random_count_range[0], random_count_range[1])
        for mz in mz_values
    ]

    plt.bar(mz_values_list, intensities, width=1)
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.show()

def fragment_peptide(peptide):
    """
    Determines the b- and y-ion sequences of a given peptide

    Parameters
    ----------
    peptide: str
        Amino acid sequence

    Returns
    -------
    dict
        Lists of the b and y-ion sequences
    """
    b_ion_list = []
    y_ion_list = []

    for i in range(1, len(peptide)):
        b_ion_list.append(peptide[:i])
        y_ion_list.insert(0, (peptide[i:]))

    return {"b_ion_list": b_ion_list, 
            "y_ion_list": y_ion_list
            }
        
