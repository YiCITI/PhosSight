import os
import sys
import re


# read fasta file
def read_fasta(file, regular):

    # check the file path
    if not os.path.exists(file):  #判断file这个文件存不存在，存在返回True，不存在返回False
        print("Error: The input file does not exist. Please check again.")
        sys.exit(1)

    # record protein names and sequences
    with open(file) as fasta:
        records = fasta.read()

    # FASTA file must start with character '>'
    if re.search('>', records) == None:   #未匹配到
        print("Error: The input file seems not in FASTA format!")
        sys.exit(1)

    # check the regular expression
    elif re.search(regular, records) == None:
        print("Error: Cannot parse the fasta file by the regular expression.")
        sys.exit(1)
    
    # fasta list
    records = re.split('(>.*?)\\n', records)[1:]   #返回一个列表，其中字符串在每次匹配时被拆分
    #['>Packet: PTM_Trainingskit_Kmod', 'PEPTIDEKIVFSGNLFQHQEDSKKLQDEIQNMKEEMARI...PEPTIDE\n',
    # '>Packet: PTM_Trainingskit_Pmod', 'PEPTIDEK.....',...,]

    length = len(records)
    fasta_list = []

    for ind in range(0, length, 2):
        name = re.split(regular, records[ind])[1]       #获取该条fasta的名称(以>开头，空格结尾)

        sequence = records[ind + 1].replace('\n', '')   #获取该条fasta的序列

        fasta_list += [[name, sequence]]

    return fasta_list