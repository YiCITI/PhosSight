import pyarrow.parquet as pq
import pyarrow.compute as pc
import pyarrow as pa # Added for pa.array and pa.string()
from tqdm import tqdm
import os
from Bio import SeqIO



def process_peptide_to_modified_sequence(peptide_sequence: str) -> tuple[str, int]:
    """
    Transforms a peptide sequence with lowercase s,t,y for phosphorylation
    into a modified sequence string (e.g., S(UniMod:21)) and counts phospho sites.
    The uppercase letter C is followed by (UniMod:4) for cysteine carbamidomethylation.

    Args:
        peptide_sequence (str): The input peptide sequence.

    Returns:
        tuple[str, int]: The transformed peptide sequence and the count of phosphorylation sites.
    """
    modified_sequence_parts = []
    phospho_count = 0
    for char_code in peptide_sequence:
        if char_code == 's':
            modified_sequence_parts.append("S(UniMod:21)")
            phospho_count += 1
        elif char_code == 't':
            modified_sequence_parts.append("T(UniMod:21)")
            phospho_count += 1
        elif char_code == 'y':
            modified_sequence_parts.append("Y(UniMod:21)")
            phospho_count += 1
        elif 'a' <= char_code <= 'z':  # Other lowercase letters
            print(f"Warning: Peptide '{peptide_sequence}' contains lowercase letter '{char_code}' which is not s, t, or y. Converting to uppercase.")
            modified_sequence_parts.append(char_code.upper())
        elif char_code == 'C':
            modified_sequence_parts.append("C(UniMod:4)")
        else:  # Uppercase letters or non-alphabetic characters
            modified_sequence_parts.append(char_code)
    return "".join(modified_sequence_parts), phospho_count


def _read_and_process_peptides_from_txt(txt_peptide_file_path: str) -> set[str]:
    """
    Reads peptides from a TXT file, processes them, and returns a set of modified sequences.

    Args:
        txt_peptide_file_path (str): Path to the input TXT file with peptide sequences.

    Returns:
        set[str]: A set of unique processed peptide sequences.
                 Returns an empty set if the file is not found or an error occurs.
    """
    processed_peptides_set = set()
    print(f"Processing peptides from: {txt_peptide_file_path}")
    try:
        with open(txt_peptide_file_path, 'r') as f:
            total_lines = sum(1 for _ in f)  # Count total lines for tqdm progress bar
            f.seek(0)
            # Wrap the enumerate(f, 1) with tqdm for a progress bar
            for line_num, line in tqdm(enumerate(f, 1), desc="Processing peptides", total=total_lines, unit="line"):
                original_peptide = line.strip()
                if not original_peptide:
                    continue

                modified_peptide, phospho_sites = process_peptide_to_modified_sequence(original_peptide)
                
                if phospho_sites > 1:
                    print(f"Warning: Peptide '{original_peptide}' (line {line_num}) has {phospho_sites} phosphorylation sites. Skipped.")
                    # And skip
                    continue

                processed_peptides_set.add(modified_peptide)

    except FileNotFoundError:
        print(f"Error: TXT peptide file not found at {txt_peptide_file_path}")
        return set()
    except Exception as e:
        print(f"Error reading or processing TXT peptide file: {e}")
        return set()
    
    return processed_peptides_set


def _read_and_process_peptides_from_fasta(fasta_peptide_file_path: str) -> set[str]:
    """
    Reads peptides from a FASTA file, processes them, and returns a set of modified sequences.

    Args:
        fasta_peptide_file_path (str): Path to the input FASTA file with peptide sequences.
        
    Returns:
        set[str]: A set of unique processed peptide sequences.
                  Returns an empty set if the file is not found or an error occurs.
    """
    processed_peptides_set = set()
    print(f"Processing peptides from FASTA: {fasta_peptide_file_path}")
    try:
        def count_fasta_records(file_path):
            with open(file_path, 'r') as f:
                return sum(1 for line in f if line.startswith('>'))
        total_records = count_fasta_records(fasta_peptide_file_path)
        
        # Use Biopython SeqIO to iterate over FASTA records
        fasta_peptide = SeqIO.parse(fasta_peptide_file_path, "fasta")
        for rec_idx, record in enumerate(tqdm(fasta_peptide, desc="Processing FASTA peptides", total=total_records), start=1):
            original_peptide = str(record.seq).strip()
            if not original_peptide:
                continue

            modified_peptide, phospho_sites = process_peptide_to_modified_sequence(original_peptide)
            if phospho_sites > 1:
                print(f"Warning: Peptide '{original_peptide}' (record {rec_idx}) has {phospho_sites} phosphorylation sites. Skipped.")
                continue

            processed_peptides_set.add(modified_peptide)

    except FileNotFoundError:
        print(f"Error: FASTA peptide file not found at {fasta_peptide_file_path}")
        return set()
    except Exception as e:
        print(f"Error reading or processing FASTA peptide file: {e}")
        return set()

    return processed_peptides_set


def _filter_parquet_by_peptide_set(
    peptide_set: set[str],
    input_parquet_path: str,
    output_parquet_path: str,
    parquet_col_to_filter: str = "Modified.Sequence"
):
    """
    Filters a Parquet file based on a set of processed peptides, saving the result.

    Args:
        peptide_set (set[str]): Set of processed peptide sequences.
        input_parquet_path (str): Path to the input Parquet spectral library.
        output_parquet_path (str): Path to save the filtered Parquet data.
        parquet_col_to_filter (str): Column name in Parquet to filter on.
    """
    if not peptide_set:
        print("Error: Peptide set is empty. Cannot proceed with Parquet filtering based on an empty list.")
        return

    print(f"Total unique processed peptide sequences for filtering: {len(peptide_set)}")

    try:
        table = pq.read_table(input_parquet_path)
        print(f"Read Parquet file: {input_parquet_path} with {len(table)} rows.")

        if parquet_col_to_filter not in table.column_names:
            raise KeyError(f"Error: Column '{parquet_col_to_filter}' not found in Parquet file. Available columns: {table.column_names}")

        sequences_in_parquet = table[parquet_col_to_filter]

        # Convert the set of peptides to a PyArrow Array for pc.is_in
        value_set_array = pa.array(list(peptide_set), type=pa.string())
        mask = pc.is_in(sequences_in_parquet, value_set=value_set_array)
        filtered_table = table.filter(mask)
        
        # Ensure output directory exists
        output_dir = os.path.dirname(output_parquet_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            print(f"Created output directory: {output_dir}")
        
        pq.write_table(filtered_table, output_parquet_path)
        print(f"Filtering complete. Output saved to: {output_parquet_path}")
        print(f"Original Parquet row count: {len(table)}")
        print(f"Filtered Parquet row count: {len(filtered_table)}")

    except FileNotFoundError:
        print(f"Error: Input Parquet file not found at {input_parquet_path}")
    except PermissionError:
        print(f"Error: Permission denied when writing to {output_parquet_path}")
    except Exception as e:
        print(f"Error processing Parquet file: {e}")


def filter_parquet_by_peptide_list(
    txt_peptide_file_path: str,
    input_parquet_path: str,
    output_parquet_path: str,
    parquet_col_to_filter: str = "Modified.Sequence"
):
    """
    Reads peptides from a TXT file, processes them for UniMod:21 (phosphorylation) 
    modifications and UniMod:4 (carbamidomethylation) modifications,
    then filters a Parquet file based on these processed peptides, saving the result.

    Args:
        txt_peptide_file_path (str): Path to the input TXT file with peptide sequences.
        input_parquet_path (str): Path to the input Parquet spectral library.
        output_parquet_path (str): Path to save the filtered Parquet data.
        parquet_col_to_filter (str): Column name in Parquet to filter on.
    """
    processed_peptides_set = _read_and_process_peptides_from_txt(txt_peptide_file_path)

    _filter_parquet_by_peptide_set(
        peptide_set=processed_peptides_set,
        input_parquet_path=input_parquet_path,
        output_parquet_path=output_parquet_path,
        parquet_col_to_filter=parquet_col_to_filter
    )
    
    
def filter_parquet_by_peptide_fasta(
    fasta_peptide_file_path: str,
    input_parquet_path: str,
    output_parquet_path: str,
    parquet_col_to_filter: str = "Modified.Sequence"
):
    """
    Reads peptides from a FASTA file, processes them for UniMod:21 (phosphorylation) 
    modifications and UniMod:4 (carbamidomethylation) modifications,
    then filters a Parquet file based on these processed peptides, saving the result.

    Args:
        fasta_peptide_file_path (str): Path to the input FASTA file with peptide sequences.
        input_parquet_path (str): Path to the input Parquet spectral library.
        output_parquet_path (str): Path to save the filtered Parquet data.
        parquet_col_to_filter (str): Column name in Parquet to filter on.
    """
    processed_peptides_set = _read_and_process_peptides_from_fasta(fasta_peptide_file_path)

    _filter_parquet_by_peptide_set(
        peptide_set=processed_peptides_set,
        input_parquet_path=input_parquet_path,
        output_parquet_path=output_parquet_path,
        parquet_col_to_filter=parquet_col_to_filter
)
    
        
def exclude_syn_pep_variants_in_parquet_by_peptide_list(
    txt_syn_pep_file_path: str,
    input_parquet_path: str,
    output_parquet_path: str,
    parquet_col_to_filter: str = "Modified.Sequence",
    protein_col_to_filter: str = "Protein.Ids"
):
    """
    Reads peptides from a TXT file, processes them for UniMod:21 (phosphorylation) 
    modifications and UniMod:4 (carbamidomethylation) modifications,
    then filters a Parquet file based on these processed peptides, saving the result.

    Args:
        txt_syn_pep_file_path (str): Path to the input TXT file with peptide sequences.
        input_parquet_path (str): Path to the input Parquet spectral library.
        output_parquet_path (str): Path to save the filtered Parquet data.
        parquet_col_to_filter (str): Column name in Parquet to filter on.
        protein_col_to_filter (str): Column name in Parquet to filter proteins (search for synthetic proteins/peptides).
    """
    processed_peptides_set = _read_and_process_peptides_from_txt(txt_syn_pep_file_path)

    if not processed_peptides_set:
        print("Error: Peptide list is empty (e.g., file not found, file empty, or error during processing). Cannot proceed with Parquet filtering based on an empty list.")
        return

    print(f"Total unique processed peptide sequences for filtering: {len(processed_peptides_set)}")

    try:
        table = pq.read_table(input_parquet_path)
        print(f"Read Parquet file: {input_parquet_path} with {len(table)} rows.")

        if parquet_col_to_filter not in table.column_names:
            raise KeyError(f"Error: Column '{parquet_col_to_filter}' not found in Parquet file. Available columns: {table.column_names}")

        sequences_in_parquet = table[parquet_col_to_filter]
        protein_ids_in_parquet = table[protein_col_to_filter]

        # Create filter conditions: 
        #   keep rows if peptide matches target list OR protein ID doesn't contain "syn"
        value_set_array = pa.array(processed_peptides_set, type=pa.string())
        peptide_match_mask = pc.is_in(sequences_in_parquet, value_set=value_set_array)
        protein_not_syn_mask = pc.invert(pc.match_substring(protein_ids_in_parquet, "syn"))
        mask = pc.or_(peptide_match_mask, protein_not_syn_mask)
        filtered_table = table.filter(mask)
        
        # Ensure output directory exists
        output_dir = os.path.dirname(output_parquet_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            print(f"Created output directory: {output_dir}")
        
        pq.write_table(filtered_table, output_parquet_path)
        print(f"Filtering complete. Output saved to: {output_parquet_path}")
        print(f"Original Parquet row count: {len(table)}")
        print(f"Filtered Parquet row count: {len(filtered_table)}")

    except FileNotFoundError:
        print(f"Error: Input Parquet file not found at {input_parquet_path}")
    except PermissionError:
        print(f"Error: Permission denied when writing to {output_parquet_path}")
    except Exception as e:
        print(f"Error processing Parquet file: {e}")

# if __name__ == "__main__":
#     import os
#     # Define file paths (replace with your actual paths)
#     script_dir = os.path.dirname(os.path.abspath(__file__))
#     test_data_dir = os.path.join(script_dir, "..", "test") # Assumes test data is in spec_parquet_filter/test/

#     txt_file = os.path.join(test_data_dir, "peptides_for_filtering.txt")
#     parquet_in = os.path.join(test_data_dir, "test.parquet") # Using your existing test.parquet
#     parquet_out = os.path.join(test_data_dir, "filtered_by_pep_list.parquet")
    
#     filter_parquet_by_peptide_list(txt_file, parquet_in, parquet_out)