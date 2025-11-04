# hyper-parameters
nearnum = 15
illegal_mer = 'BJOUX'


# get the left and right mers of the sites
def left_and_right_mer(pro_seq, ed_ind, cut_sites):
    # pro_seq    PEPTIDEKIVFSGNLFQHQEDSKKLQDEIQNMKEEMARIVSGKDYNVTANS.....
    # ed_ind     1 or 2 or 3
    # cut_sites  [-1, 7, 22, 23, 32, 37, 42, 51, 58, 63,....]
    global nearnum  #15
    
    # left mers
    left_len = len(pro_seq[ : cut_sites[ed_ind]])   #7,22,23
    if left_len >= nearnum:   #>=15
        left_mer = pro_seq[cut_sites[ed_ind] - nearnum :cut_sites[ed_ind] + 1]
        #KIVFSGNLFQHQEDSK
        #IVFSGNLFQHQEDSKK
    else:
        left_mer = 'Z' * (nearnum - left_len) +pro_seq[ : cut_sites[ed_ind] + 1]
        #ZZZZZZZZPEPRIDEK

    # right mers
    right_len = len(pro_seq[cut_sites[ed_ind] + 1 : ])  #
    if right_len >= nearnum:
        right_mer = pro_seq[cut_sites[ed_ind] + 1 :cut_sites[ed_ind] + 1 + nearnum]
        #IVFSGNLFQHQEDSK
        #KLQDEIQNMKEEMAR
        #LQDEIQNMKEEMARI
    else:
        right_mer = pro_seq[cut_sites[ed_ind] + 1 : ] +'Z' * (nearnum - right_len)
        #
    
    # left_mer = '...K/R', right_mer = '...',
    # then use full_mer to get left_mer + right_mer = '...K/R...'
    return left_mer, right_mer


# get the full 31-mer of the sites
def full_mer(pro_seq, ed_ind, cut_sites):
    # pro_seq    PEPTIDEKIVFSGNLFQHQEDSKKLQDEIQNMKEEMARIVSGKDYNVTANS.....
    # ed_ind     0 or 1 or 2 or 3
    # cut_sites  [-1, 7, 22, 23, 32, 37, 42, 51, 58, 63,....]
    if all([ed_ind > 0,ed_ind < len(cut_sites) - 1,cut_sites[ed_ind] != cut_sites[0]]):  #all([里面的元素都不为0或空])
        left_mer, right_mer = left_and_right_mer(pro_seq, ed_ind, cut_sites)
        mer = left_mer + right_mer
    else:
        mer = '*'
    return mer


# get the N-terminal digested peptides and the corresponding mers
def nterminal_pep_and_mers(pro_seq, cut_sites, terminal, missed_cleavages,min_len, max_len):
    #pro_seq    PEPTIDEKIVFSGNLFQHQEDSKKLQDEIQNMKEEMARIVSGKDYNVTANS.....
    #cut_sites  [-1, 7, 22, 23, 32, 37, 42, 51, 58, 63,....]
    #terminal   C
    #missed_cleavages 2
    #min_len    7
    #max_len    47
    nterm_seqs = []
    # print(len(cut_sites))   404

    for add_ind in range(1, len(cut_sites)):

        # check the number of missed cleavage sites
        if add_ind >= missed_cleavages + 2:  #add_ind=1,2,3
            break

        # the digested peptide with different cutting terminals
        if terminal =='C':
            pep = pro_seq[cut_sites[0] + 1 : cut_sites[add_ind] + 1]   #PEPTIDEK,PEPTIDEKIVFSGNLFQHQEDSK,PEPTIDEKIVFSGNLFQHQEDSKK
        else:
            pep = pro_seq[cut_sites[0] : cut_sites[add_ind]]

        # get the corresponding mers of the N-terminal digested peptide
        lpep = len(pep)  #8,23,24

# =============================================================================
#         If the protein N-terminal amino acid is a candidate site and cutting
#         terminal is the N-term, then the digested peptide will be '' and the
#         whole loop should be broken. Otherwise, the N-terminal digested 
#         peptides will be generated twice.
# =============================================================================
        if pep == '':
            break
# =============================================================================
#         The N-terminal digested peptides starting with M have two possible 
#         cases, e.g., MPEPTIDESK | PEPTIDESK, and the digestibility of the 
#         C-terminal cutting site is also set to 1 even if the site is 'M'.
# =============================================================================
        if lpep >= min_len and lpep <= max_len + 1:             # 48>=lpep>=7

            # 31-mer of left site
            left_mer = full_mer(pro_seq, 0, cut_sites)          #*
            if len(set(illegal_mer) - set(left_mer)) != 5:      #不进入 ==5
                continue

            # 31-mer of right site
            right_mer = full_mer(pro_seq, add_ind, cut_sites)
            # ZZZZZZZZPEPRIDEK  IVFSGNLFQHQEDSK
            # KIVFSGNLFQHQEDSK  KLQDEIQNMKEEMAR
            # IVFSGNLFQHQEDSKK  LQDEIQNMKEEMARI
            if len(set(illegal_mer) - set(right_mer)) != 5:     #不进入 ==5
                continue

            missed_mers =''
            # if add_ind > 1, then there is/are cleavage site/sites
            if add_ind > 1:  #2,3
                # 31-mers of missed cleavage sites
                # missed_num is the number of missed cleavage sites
                for missed_num in range(1, add_ind):  #1, 1 2
                    missed_left, missed_right = left_and_right_mer(pro_seq, missed_num, cut_sites)
                    missed_mer = missed_left + missed_right
                    missed_mers += (missed_mer + ',')
            #ZZZZZZZZPEPTIDEKIVFSGNLFQHQEDSK,
            #ZZZZZZZZPEPTIDEKIVFSGNLFQHQEDSK,KIVFSGNLFQHQEDSKKLQDEIQNMKEEMAR,


            # check illegal amino acids
            if len(set(illegal_mer) - set(missed_mers)) != 5:  #不进入 ==5
                continue

            mers = left_mer + '\t' + right_mer + '\t' + missed_mers.rstrip(',')
            # *       ZZZZZZZZPEPTIDEKIVFSGNLFQHQEDSK
            # *       KIVFSGNLFQHQEDSKKLQDEIQNMKEEMAR ZZZZZZZZPEPTIDEKIVFSGNLFQHQEDSK
            # *       IVFSGNLFQHQEDSKKLQDEIQNMKEEMARI ZZZZZZZZPEPTIDEKIVFSGNLFQHQEDSK,KIVFSGNLFQHQEDSKKLQDEIQNMKEEMAR

            # check the two cases starting with and without 'M'
            if lpep <= max_len:
                nterm_seqs += [pep + '\t' + mers]
            if pep[0] == 'M' and lpep - 1 >= min_len:
                nterm_seqs += [pep[1:] + '\t' + mers]
    
    # [peptide, left_mer, right_mer, missed_mers]
    # print(nterm_seqs) #['PEPTIDEK\t*\tZZZZZZZZPEPTIDEKIVFSGNLFQHQEDSK\t',
    #                   'PEPTIDEKIVFSGNLFQHQEDSK\t*\tKIVFSGNLFQHQEDSKKLQDEIQNMKEEMAR\tZZZZZZZZPEPTIDEKIVFSGNLFQHQEDSK',
    #                   'PEPTIDEKIVFSGNLFQHQEDSKK\t*\tIVFSGNLFQHQEDSKKLQDEIQNMKEEMARI\tZZZZZZZZPEPTIDEKIVFSGNLFQHQEDSK,KIVFSGNLFQHQEDSKKLQDEIQNMKEEMAR']
    return nterm_seqs


# get all theoretical peptides and corresponding 31-mers for each protein
def peps_and_mers(pro_seq, cut_sites, terminal, missed_cleavages, min_len, max_len):
    '''
    Parameters
    ----------
    pro_seq : str   PEPTIDEKIVFSGNLFQHQEDSKKLQDEIQNMKEEMARIVSGKDYNVTANS.....
        The protein sequence.
    cut_sites : list  [-1, 7, 22, 23, 32, 37, 42, 51, 58, 63,....]
        The locations of all the candidate sites in the protein sequence.
    terminal : str    C or N
        The cutting terminal.
    missed_cleavages : int   2
        The maximum number of missed cleavage sites allowed in the peptides.
    min_len : int   7
        The minimum length of digested peptides.
    max_len : int   47
        The maximum length of digested peptides.
    
    Returns
    -------
    seqs : list
        The list of digested peptides with corresponding 31-mers.
    '''

    # the N-terminal digested peptides and the corresponding mers
    seqs = nterminal_pep_and_mers(pro_seq, cut_sites, terminal, missed_cleavages, min_len, max_len)
    # ['PEPTIDEK\t*\tZZZZZZZZPEPTIDEKIVFSGNLFQHQEDSK\t',
    #  'PEPTIDEKIVFSGNLFQHQEDSK\t*\tKIVFSGNLFQHQEDSKKLQDEIQNMKEEMAR\tZZZZZZZZPEPTIDEKIVFSGNLFQHQEDSK',
    #  'PEPTIDEKIVFSGNLFQHQEDSKK\t *\tIVFSGNLFQHQEDSKKLQDEIQNMKEEMARI\tZZZZZZZZPEPTIDEKIVFSGNLFQHQEDSK, KIVFSGNLFQHQEDSKKLQDEIQNMKEEMAR']

    len_cut = len(cut_sites)  #404

    # index of the left site's position (started with the first site)
    for st_ind in range(1, len_cut - 1):  #遍历cut_sites的长度1-402

        # the number of sites on the peptide
        for add_ind in range(1, len_cut - st_ind): #1，403；1，402；1，401；..

            # check the number of missed cleavage sites
            if add_ind >= missed_cleavages + 2:   #>=4(1,2,3)
                break

            # the digested peptide with different cutting terminals
            if terminal =='C':
                pep = pro_seq[cut_sites[st_ind] + 1 :cut_sites[st_ind + add_ind] + 1]
            else:
                pep = pro_seq[cut_sites[st_ind] :cut_sites[st_ind + add_ind]]
            
            # get the corresponding mers of the digested peptide
            lpep = len(pep)
            if lpep >= min_len and lpep <= max_len:
                # 31-mer of left site
                left_mer = full_mer(pro_seq, st_ind, cut_sites)
                if len(set(illegal_mer) - set(left_mer)) != 5:
                    continue
                # 31-mer of right site
                right_mer = full_mer(pro_seq, st_ind + add_ind, cut_sites)
                if len(set(illegal_mer) - set(right_mer)) != 5:
                    continue
                missed_mers = ''
                # if add_ind > 1, then there is/are cleavage site/sites
                if add_ind > 1:
                    # 31-mers of missed cleavage sites
                    for missed_num in range(1, add_ind):
                    # missed_num is the number of missed cleavage sites
                        missed_left, missed_right =left_and_right_mer(pro_seq, st_ind + missed_num,cut_sites)
                        missed_mer = missed_left + missed_right
                        missed_mers += (missed_mer + ',')
                # check illegal amino acids
                if len(set(illegal_mer) - set(missed_mers)) != 5:
                    continue
                seqs += [pep + '\t' + left_mer + '\t' + right_mer +'\t' + missed_mers.rstrip(',')]
    
    # [peptide, left_mer, right_mer, missed_mers]
    # print(seqs)
    return seqs


# protein --> peptide | left 31-mer | right 31-mer | missed 31-mers
def digestion(pro_seq, sites, terminal, missed_cleavages, min_len, max_len):
    #pro_seq是fasta里面的一条蛋白质序列
    #PEPTIDEKIVFSGNLFQHQEDSKKLQDEIQNMKEEMARIVSGKDYNVTANS.....
    lpro = len(pro_seq)

    # indexes of the candidate cleavage sites
    cut_sites = [ind for ind in range(lpro) if pro_seq[ind] in sites]  #从0开始
    #[7, 22, 23, 32, 37, 42, 51, 58, 63, 68,...]

    # indexes of the start and end positions
# =============================================================================
#     There are some differences between C-terminal and N-terminal cleavage, 
#     but both need to start at the first amino acid and end at the last one.
# ============================================================================+
    if terminal == 'C':
        cut_sites.insert(0, -1)        #在cut_sites开头插入-1
        if cut_sites[-1] != lpro - 1:  #如果最后一个切割位点不等于该fasta长度减1,则将最后一个也算为切割位点
            cut_sites += [lpro - 1]
    else:
        cut_sites.insert(0, 0)         #在cut_sites开头插入0
        cut_sites += [lpro]

    # check if there is any candidate sites to cut
    digested_seqs = []
    if len(cut_sites) > 2:
        digested_seqs += peps_and_mers(pro_seq, cut_sites, terminal, missed_cleavages, min_len, max_len)

    elif lpro >= min_len and lpro <= max_len + 1:
        print("Warning: Protein sequence %s has no candidate sites!" % pro_seq)
        if lpro <= max_len:
            digested_seqs += [pro_seq]
        if pro_seq[0] == 'M' and lpro - 1 >= min_len:
            digested_seqs += [pro_seq[1:]]
    
    return digested_seqs
