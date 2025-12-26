import os
from typing import Dict
from Bio import SeqIO
import pandas as pd
import logging

logger = logging.getLogger(__name__)

class ProteinSpeciesMapper:
    """
    负责从多个fasta文件加载蛋白质ID与物种的映射，并提供ID到物种的查询接口。
    """
    def __init__(self, species_fasta_dict: Dict[str, str]):
        """
        species_fasta_dict: 物种名称(str) -> fasta文件路径(str) 的字典
        """
        self.species_fasta_dict = species_fasta_dict
        self.protein_species_map = self._load_protein_species_map()
        logger.info(f"已加载 {len(self.protein_species_map)} 个蛋白质ID的物种映射关系")


    def _load_protein_species_map(self) -> Dict[str, str]:
        protein_species_map = {}
        for species_name, fasta_file in self.species_fasta_dict.items():
            species_map = self._get_protein_one_species(fasta_file, species_name)
            protein_species_map.update(species_map)
        return protein_species_map


    def _get_protein_one_species(self, fasta_file: str, species_name: str) -> Dict[str, str]:
        protein_species_map = {}
        if not os.path.exists(fasta_file):
            logger.warning(f"Fasta文件不存在: {fasta_file}")
            return protein_species_map
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                parts = record.id.split('|')
                if len(parts) >= 3:
                    accession = parts[1]
                else:
                    accession = record.id
                protein_species_map[accession] = species_name
        logger.info(f"从 {fasta_file} 加载了 {len(protein_species_map)} 个 {species_name} 蛋白质")
        return protein_species_map


    def map_ids_to_species(self, protein_ids_str: str) -> str:
        """
        将蛋白质ID字符串映射为物种，处理多个ID的情况
        """
        if pd.isna(protein_ids_str):
            raise ValueError("蛋白质ID为空")
        protein_ids = [pid.strip() for pid in str(protein_ids_str).split(';')]
        species_set = set()
        missing_ids = []
        for protein_id in protein_ids:
            if protein_id in self.protein_species_map:
                species_set.add(self.protein_species_map[protein_id])
            else:
                missing_ids.append(protein_id)
        if missing_ids:
            raise ValueError(f"以下蛋白质ID在映射字典中未找到: {missing_ids}")
        if len(species_set) > 1:
            # 多物种时用分号连接所有物种名称
            return ';'.join(sorted(species_set))
        elif len(species_set) == 0:
            raise ValueError(f"蛋白质ID '{protein_ids_str}' 未找到对应物种")
        return species_set.pop()
