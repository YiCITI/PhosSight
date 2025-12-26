from process_fasta.in_silico_digestion import digestion
from Bio import SeqIO
from tqdm import tqdm

def read_fasta(file_path):
    """Read sequences from a FASTA file."""
    return [(str(record.id), str(record.seq)) for record in SeqIO.parse(file_path, "fasta")]

def in_silico_digest(sequence, miss_cleavage=4, min_length=7, max_length=46):
    """Perform in silico digestion of a sequence."""
    # Assuming digestion function takes a sequence and returns a list of peptides
    return digestion(sequence, 'KR', 'C', miss_cleavage, min_length, max_length)

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

def digest_fasta(file_path, miss_cleavage=4, min_length=7, max_length=46):
    """Digest sequences from a FASTA file and extend with phosphorylation modifications.
    Returns a list of (id, peptide_seq) tuples."""
    proteins = read_fasta(file_path)  # [(id, seq), ...]
    all_extended_peptides = []
    for id, seq in tqdm(proteins, desc="Processing sequences"):
        digest_result = in_silico_digest(seq, miss_cleavage, min_length, max_length)
        extended_peptides = extend(digest_result)
        # 绑定id和peptide
        for pep in extended_peptides:
            pep = pep.split('\t')[0] # 这里加一句
            all_extended_peptides.append((id, pep))
    return all_extended_peptides

def generate_modified_pep_fasta_file(
    file_path, 
    output_file, 
    miss_cleavage=4, 
    min_length=7, 
    max_length=46,
    file_type='Uniprot',
    is_replace_STY_to_BJX=True
):
    """Generate a new FASTA file with the digested and extended peptides, keeping id info."""
    print(f"digesting {file_path} with miss_cleavage={miss_cleavage}, min_length={min_length}, max_length={max_length}")
    digest_result = digest_fasta(file_path, miss_cleavage, min_length, max_length)  # [(id, seq), ...]
    # Only keep the first part of the sequence before any tab character
    peptides = {}
    for id, seq in tqdm(digest_result, desc="Processing peptides"):
        seq = seq.split('\t')[0]
        if seq.strip():
            if file_type == 'Uniprot':
                id = id.split('|')[1]
            elif file_type == 'IPI':
                idx_vert_line = id.find('|')
                if idx_vert_line != -1:
                    id = id[4:idx_vert_line]
                else:
                    id = id[4:]
            else:
                raise ValueError(f"Unsupported file type: {file_type}. Supported types are 'Uniprot' and 'IPI'.")
            if seq not in peptides:
                peptides[seq] = set()
            peptides[seq].add(id)

    with open(output_file, 'w') as f:
        for index, (seq, ids) in enumerate(tqdm(peptides.items(), desc="Writing peptides to FASTA")):
            
            # Skip sequences that contain 'BJX' originally to avoid confusion
            # This is to ensure that we do not process sequences that are already in BJX format
            if 'BJX' in seq:
                print(f"Skipping sequence {seq} as it contains 'BJX' which is not expected in the original sequences.")
                continue
            
            # Replace 's', 't', 'y' with 'B', 'J', 'X' respectively
            if is_replace_STY_to_BJX:
                seq = seq.replace('s', 'B').replace('t', 'J').replace('y', 'X')
            f.write(f">{'|'.join(sorted(ids))}_{index}\n{seq}\n")