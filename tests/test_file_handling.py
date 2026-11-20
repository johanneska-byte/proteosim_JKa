def test_fasta_reader(filepath):
    tmp_fasta_path = '../data/sample_proteins.fasta'
    #tmp_fasta_path = filepath

    protein_map = read_fasta(tmp_fasta_path)
    
    # Replace the strings with your fasta content
    # which you expect to be now available as a dictionary
    assert protein_map["P11802"].startswith("MATSRYEPVAEIGVGAYGTV")
    assert protein_map["A0A087WTH1"].startswith("MEDEEKAVEILGNTEAAHPP")
