import sys

def txt_to_fasta(txt_path, fasta_path, prefix):
    """Convert a text file of sequences to a FASTA file.

    Parameters
    ----------
    txt_path : str
        Path to the input text file containing sequences, one per line.
    fasta_path : str
        Path to the output FASTA file.
    prefix : str
        Prefix to use for the protein names in the FASTA file.
    """
    with open(txt_path, 'r', encoding='utf-8') as txt_file, open(fasta_path, 'w', encoding='utf-8') as fasta_file:
        for idx, line in enumerate(txt_file):
            seq = line.strip()
            if seq:
                protein_name = f"{prefix}_{idx}"
                fasta_file.write(f">{protein_name}\n{seq}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python txt2fasta.py <input_txt> <output_fasta> <prefix>")
        sys.exit(1)
    txt_to_fasta(sys.argv[1], sys.argv[2], sys.argv[3])