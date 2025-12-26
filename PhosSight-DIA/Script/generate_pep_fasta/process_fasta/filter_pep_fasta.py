from Bio import SeqIO
from os import path
import re

def filter_peptides(filtered_peptides_list_file, all_peptides_fasta_file, filtered_peptides_fasta_file):
    """
    Filters peptides from the input file and writes them to the output file.
    
    Parameters:
    filtered_peptides_list_file (str): Path to the text file containing peptides to keep.
    all_peptides_fasta_file (str): Path to the input FASTA file containing all peptides.
    filtered_peptides_fasta_file (str): Path to the output FASTA file for filtered peptides.
    """
    # Read the peptides for filtering from the text file
    with open(filtered_peptides_list_file) as f:
        keep_set = set(line.strip() for line in f if line.strip())

    # Filter the FASTA file
    with open(all_peptides_fasta_file) as fasta_in, open(filtered_peptides_fasta_file, "w") as fasta_out:
        for record in SeqIO.parse(fasta_in, "fasta"):
            if str(record.seq) in keep_set:
                SeqIO.write(record, fasta_out, "fasta")
                


def filter_peptides_above_thres(
    peptides_list_file_for_filtering, 
    all_peptides_fasta_file,
    output_txt_file, 
    score_threshold
    ):
    """
    Filters peptides from the input file based on a score threshold and writes them to the output file.
    
    Parameters:
    peptides_list_file_for_filtering (str): Path to the text file containing peptides and their scores.
    all_peptides_fasta_file (str): Path to the input FASTA file containing all peptides.
    output_txt_file (str): Path to the output text file for filtered peptides.
    score_threshold (float): Minimum score threshold for filtering peptides.
    """
    # Read the peptides for filtering from the text file
    keep_set = set()
    with open(peptides_list_file_for_filtering) as f:
        for line in f:
            parts = re.split(r'[,\t]', line.strip())
            if len(parts) == 2:
                peptide, score_str = parts
                try:
                    score = float(score_str)
                    if score >= score_threshold:
                        keep_set.add(peptide)
                except ValueError:
                    continue

    # Filter the FASTA file
    with open(all_peptides_fasta_file) as fasta_in, open(output_txt_file, "w") as txt_out:
        for record in SeqIO.parse(fasta_in, "fasta"):
            if str(record.seq) in keep_set:
                txt_out.write(f"{record.seq}\n")


def filter_peptides_by_ratio(
    peptides_list_file_for_filtering, 
    all_peptides_fasta_file,
    output_txt_file, 
    ratio_to_keep
    ):
    """
    Sort peptides by score and keep the top ratio_to_keep fraction, then write to output txt.

    Parameters:
    peptides_list_file_for_filtering (str): Path to the text file containing peptides and scores.
    all_peptides_fasta_file (str): Path to the input FASTA file.
    output_txt_file (str): Path to the output text file.
    ratio_to_keep (float): Fraction to keep (between 0 and 1).
    """
    peptide_score_list = []
    with open(peptides_list_file_for_filtering) as f:
        for line in f:
            parts = re.split(r'[,\t]', line.strip())
            if len(parts) == 2:
                peptide, score_str = parts
                try:
                    score = float(score_str)
                    peptide_score_list.append((peptide, score))
                except ValueError:
                    continue

    peptide_score_list.sort(key=lambda x: x[1], reverse=True)
    total_peptides = len(peptide_score_list)
    num_to_keep = int(total_peptides * ratio_to_keep)
    keep_set = set(peptide for peptide, score in peptide_score_list[:num_to_keep])

    with open(all_peptides_fasta_file) as fasta_in, open(output_txt_file, "w") as txt_out:
        for record in SeqIO.parse(fasta_in, "fasta"):
            if str(record.seq) in keep_set:
                txt_out.write(f"{record.seq}\n")


if __name__ == "__main__":
    
    base_dir = path.dirname(path.abspath(__file__))
    
    filtered_peptides_list_file = path.join(base_dir, "filtered_peptides_BJX.txt")
    all_peptides_fasta_file = path.join(base_dir, "all_peptides.fasta")
    filtered_peptides_fasta_file = path.join(base_dir, "filtered_pep_pretrained.fasta")

