def test_digest_protein_sequence():
    dummy_seq_1 = "MTSRSLANLQRRRQQQQQQRAFSLVGKGGFGKVYKATNKLTDDHVNIKIVDFGSCAKLSLVGTPNWMAPEVIRGKAYGPEVDIWSLGVILYEMAAGYPPFFADQPIQIYQVSRGEGLLGVVNYLHKGSLTFWQLEQKRGPNQPNLKSPQAQQAQ"
    dummy_seq_2 = "MDSRMASHLSSIPSLPPKLLLLGLCVGAGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHYYR"
    
    assert digest_protein_sequence(dummy_seq_1, '(?<=R)') == ['SLANLQR', 'QQQQQQR', 'GEGLLGVVNYLHKGSLTFWQLEQKR', 'GPNQPNLKSPQAQQAQ']
    assert digest_protein_sequence(dummy_seq_2, '(?<=R)') == ['KQVVIDGETCLLDILDTAGQEEYSAMR', 'DQYMR', 'TGEGFLCVFAINNTKSFEDIHYYR']
    assert digest_protein_sequence(dummy_seq_1, '(?<=[KR])(?!P)') == ['SLANLQR', 'QQQQQQR', 'AFSLVGK', 'GGFGK', 'LTDDHVNIK', 'IVDFGSCAK', 'LSLVGTPNWMAPEVIR', 'GEGLLGVVNYLHK', 'GSLTFWQLEQK', 'GPNQPNLK', 'SPQAQQAQ']
    assert digest_protein_sequence(dummy_seq_2, '(?<=[KR])(?!P)') == ['MASHLSSIPSLPPK', 'LLLLGLCVGAGVGK', 'SALTIQLIQNHFVDEYDPTIEDSYR', 'QVVIDGETCLLDILDTAGQEEYSAMR', 'DQYMR', 'TGEGFLCVFAINNTK', 'SFEDIHYYR']
    assert digest_protein_sequence(dummy_seq_2, '(?<=[KR])(?!P)', min_pep_len = 3, max_pep_len = 15) == ['MDSR', 'MASHLSSIPSLPPK', 'LLLLGLCVGAGVGK', 'DQYMR', 'TGEGFLCVFAINNTK', 'SFEDIHYYR']


def test_digest_protein_collection():

    # Cleave after aromatic amino acids (F, W, Y)
    cleave_pattern_aromatic = r'(?<=[FWY])'
    # Cleave before acidic amino acids (D, E)
    cleave_pattern_acidic = r'(?=[DE])'
    # Cleave after two consecutive glutamines (QQ)
    cleave_pattern_double_q = r'(?<=QQ)'
    
    dummy_proteins = {
    "Protein_1": "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE",
    "Protein_2": "MSPQTETKASVGFKAGVKDYKLTYYTPDYETKDTDILAAFR",
    "Protein_3": "MADQLTEEQIAEFKEAFSLFDKDGDGTITTKEGVVVR",
    "Protein_4": "MNNNKQQQQQQRGEGLLGVVNYLHKGSLTFWQLEQKRGPNQPNLKSPQAQQAQ",
    "Protein_5": "MPEPTIDERKAAAAAKPEPTIDERKLLLLLRPEPTIDE"
    }

    assert(digest_protein_collection(dummy_proteins, cleave_pattern_aromatic)) == {'Protein_1': ['ISLLF', 'SRGVF', 'RRDTHKSEIAHRF', 'KDLGE'], 'Protein_2': ['MSPQTETKASVGF', 'KAGVKDY', 'ETKDTDILAAF'], 'Protein_3': ['MADQLTEEQIAEF', 'DKDGDGTITTKEGVVVR'], 'Protein_4': ['MNNNKQQQQQQRGEGLLGVVNY', 'LHKGSLTF', 'QLEQKRGPNQPNLKSPQAQQAQ'], 'Protein_5': []}
    assert(digest_protein_collection(dummy_proteins, cleave_pattern_acidic)) == {'Protein_1': ['MKWVTFISLLFLFSSAYSRGVFRR', 'DTHKS', 'EIAHRFK'], 'Protein_2': ['MSPQT', 'ETKASVGFKAGVK', 'DYKLTYYTP', 'DILAAFR'], 'Protein_3': ['EAFSLF', 'DGTITTK', 'EGVVVR'], 'Protein_4': ['MNNNKQQQQQQRG', 'EGLLGVVNYLHKGSLTFWQL', 'EQKRGPNQPNLKSPQAQQAQ'], 'Protein_5': ['ERKAAAAAKP', 'ERKLLLLLRP']}
    assert(digest_protein_collection(dummy_proteins, cleave_pattern_double_q)) == {'Protein_1': [], 'Protein_2': [], 'Protein_3': [], 'Protein_4': ['MNNNKQQ'], 'Protein_5': []}


def test_compute_sequence_coverage():
    dummy_prot_seq = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"
    dummy_peps = [
    "MKWVTFISLLFLFSSAYSR",
    "SEIAHR",
    "FKDLGE"
    ]

    coverage = compute_sequence_coverage(dummy_prot_seq, dummy_peps)
    assert coverage == 77.5