import os
from typing import Dict
from Bio import SeqIO
import pandas as pd
import logging

logger = logging.getLogger(__name__)

class ProteinSpeciesMapper:
    """
    Load protein ID to species mapping from multiple fasta files and providing an interface for querying species by ID.
    """
    def __init__(self, species_fasta_dict: Dict[str, str]):
        """
        species_fasta_dict: Dictionary mapping species name (str) -> fasta file path (str)
        """
        self.species_fasta_dict = species_fasta_dict
        self.protein_species_map = self._load_protein_species_map()
        logger.info(f"Loaded species mapping for {len(self.protein_species_map)} protein IDs")


    def _load_protein_species_map(self) -> Dict[str, str]:
        protein_species_map = {}
        for species_name, fasta_file in self.species_fasta_dict.items():
            species_map = self._get_protein_one_species(fasta_file, species_name)
            protein_species_map.update(species_map)
        return protein_species_map


    def _get_protein_one_species(self, fasta_file: str, species_name: str) -> Dict[str, str]:
        protein_species_map = {}
        if not os.path.exists(fasta_file):
            logger.warning(f"Fasta file does not exist: {fasta_file}")
            return protein_species_map
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                parts = record.id.split('|')
                if len(parts) >= 3:
                    accession = parts[1]
                else:
                    accession = record.id
                protein_species_map[accession] = species_name
        logger.info(f"Loaded {len(protein_species_map)} {species_name} proteins from {fasta_file}")
        return protein_species_map


    def map_ids_to_species(self, protein_ids_str: str) -> str:
        """
        Map protein ID string to species, handling multiple IDs
        """
        if pd.isna(protein_ids_str):
            raise ValueError("Protein ID is empty")
        protein_ids = [pid.strip() for pid in str(protein_ids_str).split(';')]
        species_set = set()
        missing_ids = []
        for protein_id in protein_ids:
            if protein_id in self.protein_species_map:
                species_set.add(self.protein_species_map[protein_id])
            else:
                missing_ids.append(protein_id)
        if missing_ids:
            raise ValueError(f"The following protein IDs were not found in the mapping dictionary: {missing_ids}")
        if len(species_set) > 1:
            # Join all species names with semicolon for multi-species
            return ';'.join(sorted(species_set))
        elif len(species_set) == 0:
            raise ValueError(f"No corresponding species found for protein ID '{protein_ids_str}'")
        return species_set.pop()
