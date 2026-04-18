import argparse
import pandas as pd
from pathlib import Path
import re

def process_peptide(pep_str):
    if pd.isna(pep_str):
        raise ValueError("Peptide sequence cannot be NaN or missing.")
    
    pep_str = str(pep_str)
    
    mod_mapping = {
        '1': ('M', 'Oxidation of M', '15.994915'),
        '2': ('S', 'Phospho_no_loss of S', '79.966331'),
        '3': ('T', 'Phospho_no_loss of T', '79.966331'),
        '4': ('Y', 'Phospho_no_loss of Y', '79.966331')
    }
    
    new_pep_chars = []
    mods = []
    
    for char in pep_str:
        if char in mod_mapping:
            aa, mod_name, mass = mod_mapping[char]
            new_pep_chars.append(aa)
            mods.append(f"{mod_name}@{len(new_pep_chars)}[{mass}]")
        else:
            new_pep_chars.append(char)
            
    return "".join(new_pep_chars), ";".join(mods)

def main():
    parser = argparse.ArgumentParser(description="Generate PDV table and index from CSV.")
    parser.add_argument("input_csv", help="Path to the input CSV file.")
    args = parser.parse_args()
    
    input_path = Path(args.input_csv)
    
    df = pd.read_csv(input_path)
    
    # PDV table records [{"spectrum_title": ..., "peptide": ..., "charge": ..., "modification": ...}, ...]
    table_records = []
    
    # PDV index set to store unique spectrum titles
    unique_titles = set()
    
    for _, row in df.iterrows():
        spectrum_title = str(row['spectrum_title'])
        unique_titles.add(spectrum_title)
        
        charge = spectrum_title.split('.')[-1]
        
        for pep_col in ['peptide_a', 'peptide_b']:
            if pep_col in row:
                pep_seq, modifications = process_peptide(row[pep_col])
                if pep_seq:
                    table_records.append({
                        'spectrum_title': spectrum_title,
                        'peptide': pep_seq,
                        'charge': charge,
                        'modification': modifications
                    })
                else:
                    raise ValueError(f"Processed peptide sequence is empty for row with spectrum_title: {spectrum_title}")
                    
    table_df = pd.DataFrame(table_records)
    
    table_out = input_path.with_name(f"{input_path.stem}_PDV_table.txt")
    table_df.to_csv(table_out, sep='\t', index=False)
    
    index_out = input_path.with_name(f"{input_path.stem}_PDV_index.txt")
    with open(index_out, 'w') as f:
        for title in sorted(list(unique_titles)):
            f.write(f"{title}\n")
            
if __name__ == '__main__':
    main()
