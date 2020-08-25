#######检查一段核苷酸序列是否是有效序列####

Nucleotides = ['A','T','C','G']
def validataSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

#####生成随机的DNA序列#########
import random
Nucleotides = ['A','T','C','G']
randDNAStr = ''.join([random.choice(Nucleotides)
                      for nuc in range(20)])
# =============================================================================
# print(randDNAStr)
# =============================================================================
# =============================================================================
# print(validataSeq(randDNAStr))
# =============================================================================
dnaStr = validataSeq(randDNAStr)
######计算核苷酸序列中ATCG的个数#######
import collections
def countNucFrequency(seq):
# =============================================================================
#     tmpFreDict = {"A":0,"T":0,"C":0,"G":0}
#     for nuc in seq:
#         tmpFreDict[nuc] = tmpFreDict[nuc]+1
#     return tmpFreDict
# =============================================================================
    return dict(collections.Counter(seq))

# =============================================================================
# print(countNucFrequency(dnaStr))
# =============================================================================
###########transciption dna-rna##########
def transciption(seq):
    return seq.replace("T","U")
# =============================================================================
# print(transciption(dnaStr))
# =============================================================================



###################dna reverse##############
DNA_ReverseComplement = {"A":"T","T":"A","C":"G","G":"C"}
def reverse_complement(seq):
    return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]

# =============================================================================
# print(reverse_complement(dnaStr))
# =============================================================================



###########color###############
def colored(seq):
    bcolors={
        'A': '\033[92m',
        'T': '\033[91m', 
        'C': '\033[94m',
        'G': '\033[93m',
        'U': '\033[91m',
        'reset': '\033[0;0m'}
    tmpStr = ""
    for nuc in seq:
        if nuc in bcolors:
            tmpStr = tmpStr + bcolors[nuc] + nuc
        else:
            tmpStr = tmpStr + bcolors['reset']+nuc
    return tmpStr + '\033[0;0m'
print(f'\nSequence:{colored(dnaStr)}\n')
print(f'[1]+Sequence Length:{len(dnaStr)}\n')
print(colored(f'[2]+Nucleotide Frequency:{countNucFrequency(dnaStr)}\n'))
print(f'[3]+dna\rna Transcription:{colored(transciption(dnaStr))}\n')
print(f"[4]+DNA string +reverse complement:\n5'{colored(dnaStr)}3'")
print(f"  {''.join(['|' for c in range(len(dnaStr))])}")
print(f"3'{colored(reverse_complement(dnaStr))}5'\n")    
    
    
    
    
    
    
    
    