def read_seq(sequence):
    bases = {'A', 'T', 'C', 'G'}
    sequence = sequence.upper()
    
    for base in sequence:
        if base not in bases:
            print(f"Erreur : La séquence contient une base invalide '{base}'.")
            return None  

    count_for_base = {
        'A': sequence.count('A'),
        'C': sequence.count('C'),
        'G': sequence.count('G'),
        'T': sequence.count('T')
    }

    return count_for_base
