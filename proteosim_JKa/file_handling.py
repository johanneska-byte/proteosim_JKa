def read_fasta(filepath):
    """
    Reads a .fasta file, extracts Protein_IDs and connects them with associated sequence as string. Saves all Protein_ID / sequence pairs in a dictionary.

    Parameters
    ----------
    filepath:   str
                path to PASTA file
    Returns
    -------
    dictionary
        dictionary mapping ProteinIDs to corresponding sequence
    """
    
    protein_map = {}                                                                        # erzeugt leeres dictionairy
    current_id = None                                                                       # erzeugt variable 
    current_sequence = []                                                                   # erzeugt leere liste

    with open(filepath, 'r') as fasta_handle:                                             # öffnet die file unter fasta_path in einer variable namens fasta_handle
        for line in fasta_handle:                                                           # interiert über alle lines
            stripped = line.strip()                                                         # entfernt Leerzeichen
            if not stripped:                                                                # wenn Zeile nach strippen leer ist
                continue
            if stripped.startswith('>'):                                                    # wenn die Line mit > beginnt 
                if current_id is not None:                                                  # und es schon eine zugeordnete ID gibt
                    protein_map[current_id] = ''.join(current_sequence)                     # Listenelement werden zusammengefügt; dann geht das als Eintrag ins dictionairy
                    current_sequence = []                                                   # current sequence leere liste gleichsetzten
                current_id = stripped.split('|')[1]                                         # splittet die Zeichenkette bei jedem |; speichert zweites Element in Liste 
            else:
                current_sequence.append(stripped)                                           # die sequenzzeile wird in Liste eingetragen

    if current_id is not None:
        protein_map[current_id] = ''.join(current_sequence)
                                         # letztes Protein speichern; weil ja kein speichern auslösender Header mehr folgt
    n_proteins = len(protein_map)
    for prot_id, seq in protein_map.items():
        print(f'{prot_id}: {len(seq)}')

    return protein_map
