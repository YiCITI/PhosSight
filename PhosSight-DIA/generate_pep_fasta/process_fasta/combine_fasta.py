from Bio import SeqIO
from tqdm import tqdm

def combine_fasta(fasta1_path, fasta2_path, output_path):
    def get_base_name(record_id):
        # Split by last underscore to get base name
        if '_' in record_id:
            base = '_'.join(record_id.split('_')[:-1])
        else:
            base = record_id
        return base

    # Parse both fasta files
    records1 = list(SeqIO.parse(fasta1_path, "fasta"))
    records2 = list(SeqIO.parse(fasta2_path, "fasta"))

    # Combine and assign new suffixes
    combined = records1 + records2
    with open(output_path, 'w') as out:
        for idx, record in tqdm(enumerate(combined), desc="Combining FASTA files", total=len(combined)):
            base = get_base_name(record.id)
            new_id = f'{base}_{idx}'
            out.write(f'>{new_id}\n')
            out.write(f'{record.seq}\n')

# Example usage:
# combine_fasta('fasta1.fa', 'fasta2.fa', 'combined.fa')