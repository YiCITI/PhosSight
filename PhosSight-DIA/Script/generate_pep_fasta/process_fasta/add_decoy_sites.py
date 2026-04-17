import argparse
from Bio import SeqIO


def add_decoy_sites_to_fasta(input_fasta, output_fasta, prefix):
	"""Create single-site lowercase STY decoy entries from a FASTA file.

	For each input sequence, this function enumerates all positions containing
	uppercase S/T/Y and generates a new sequence where exactly one such residue
	is converted to lowercase. New entries are appended only when the sequence
	has not appeared before in the input database or in already generated decoys.
	"""
	records = list(SeqIO.parse(input_fasta, "fasta"))

	existing_sequences = set()
	for record in records:
		existing_sequences.add(str(record.seq))

	# Always enumerate from uppercase parent sequences so each generated decoy
	# contains at most one lowercase STY site.
	base_sequences = {seq.upper() for seq in existing_sequences}

	new_entries = []
	decoy_index = 0 

	for seq in base_sequences:
		for pos, aa in enumerate(seq):
			if aa in {"S", "T", "Y"}:
				decoy_seq = seq[:pos] + aa.lower() + seq[pos + 1 :]
				if decoy_seq in existing_sequences:
					continue

				decoy_id = f"{prefix}_{decoy_index}"
				new_entries.append((decoy_id, decoy_seq))
				existing_sequences.add(decoy_seq)
				decoy_index += 1

	with open(output_fasta, "w", encoding="utf-8") as out_f:
		for record in records:
			out_f.write(f">{record.id}\n{str(record.seq)}\n")
		for decoy_id, decoy_seq in new_entries:
			out_f.write(f">{decoy_id}\n{decoy_seq}\n")

	return len(new_entries)


def main():
	parser = argparse.ArgumentParser(
		description="Add single-site lowercase STY decoy peptides into FASTA"
	)
	parser.add_argument("input_fasta", type=str, help="Input FASTA path")
	parser.add_argument("output_fasta", type=str, help="Output FASTA path")
	parser.add_argument("prefix", type=str, help="Decoy sequence ID prefix")
	args = parser.parse_args()

	n_added = add_decoy_sites_to_fasta(
		input_fasta=args.input_fasta,
		output_fasta=args.output_fasta,
		prefix=args.prefix,
	)
	print(f"Added {n_added} decoy-site sequences to {args.output_fasta}")


if __name__ == "__main__":
	main()
