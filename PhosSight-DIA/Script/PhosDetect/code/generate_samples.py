from in_silico_digestion import digestion
from Bio import SeqIO

def read_fasta(file_path):
    """Read sequences from a FASTA file."""
    return [str(record.seq) for record in SeqIO.parse(file_path, "fasta")]

def in_silico_digest(sequence):
    """Perform in silico digestion of a sequence."""
    # Assuming digestion function takes a sequence and returns a list of peptides
    return digestion(sequence, 'KR', 'C', 1, 7, 30)

def extend(digest_result):
    """Make in silico phosphorylation modifications with only one modification per peptide."""
    extended_peptides = []
    for peptide in digest_result:
        extended_peptides.append(peptide)
        for i, c in enumerate(peptide):
            if c == 'S' or c == 'T' or c == 'Y':
                # change S, T or Y to s, t or y to indicate phosphorylation
                modified_peptide = peptide[:i] + c.lower() + peptide[i+1:]
                extended_peptides.append(modified_peptide)
    return list(set(extended_peptides))

def digest_fasta(file_path):
    """Digest sequences from a FASTA file and extend with phosphorylation modifications."""
    sequences = read_fasta(file_path)
    all_extended_peptides = []
    
    for sequence in sequences:
        digest_result = in_silico_digest(sequence)
        extended_peptides = extend(digest_result)
        all_extended_peptides.extend(extended_peptides)
    
    return list(set(all_extended_peptides))

