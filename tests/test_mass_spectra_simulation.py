import proteosim_JKa as ps

def test_calculate_mol_mass():
    aa_mass_dict = ps.amino_acid_mass_dalton
    test_peptide = "IINFRMATIK"
    expected_mass = 1188.50
    calculated_mass = list(ps.calculate_mol_mass(test_peptide, aa_mass_dict).values())[0]
    assert round(calculated_mass, 2) == round(expected_mass, 2)

def test_calculate_mol_mass_collection():
    aa_mass_dict = ps.amino_acid_mass_dalton
    peptides = ["AA", "PEP", "ACDK", "MKWVTF", "RLGY", "VVVV"]
    test_mol_mass_dict = {}

    for i in peptides:
        mol_mass = 0
        for amino_acid in i:
            mol_mass += ps.amino_acid_mass_dalton[amino_acid]

        test_mol_mass_dict[i] = mol_mass

    expected = test_mol_mass_dict
    actual = ps.calculate_mol_mass_collection(peptides, aa_mass_dict)

    assert actual == expected
    assert isinstance(actual, dict)
    assert len(actual) == len(peptides)

def test_calculate_mz_collection():
    peptide_mass_map = {
    "AA": 142.16,
    "PEP": 341.37,
    "ACDK": 417.49
    }
    actual_charge2 = ps.calculate_mz_collection(peptide_mass_map, charge=2)
    expected_charge2 = {'AA': {'mz': 72.087, 'mass': 142.16}, 'PEP': {'mz': 171.692, 'mass': 341.37}, 'ACDK': {'mz': 209.752, 'mass': 417.49}}
    actual_charge5 = ps.calculate_mz_collection(peptide_mass_map, charge= 5)
    expected_charge5 = {'AA': {'mz': 29.439, 'mass': 142.16}, 'PEP': {'mz': 69.281, 'mass': 341.37}, 'ACDK': {'mz': 84.50500000000001, 'mass': 417.49}}
    actual_charge2_proton_mass1 = ps.calculate_mz_collection(peptide_mass_map, charge=2, proton_mass= 1)
    expected_charge2_proton_mass1 = {'AA': {'mz': 72.08, 'mass': 142.16}, 'PEP': {'mz': 171.685, 'mass': 341.37}, 'ACDK': {'mz': 209.745, 'mass': 417.49}}

    assert actual_charge2 == expected_charge2
    assert actual_charge5 == expected_charge5
    assert actual_charge2_proton_mass1 == expected_charge2_proton_mass1

def test_fragment_peptide():
    peptide = "PEPTIDE"
    expected = {'b_ion_list': ['P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID'],
                'y_ion_list': ['EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']
                }
    actual = ps.fragment_peptide("PEPTIDE")

    assert set(actual) == set(expected)
    assert 'PEPTIDE' not in expected  # which string should not be in the expected?
