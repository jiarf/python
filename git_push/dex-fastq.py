#!/home/devdata/Software/soft/anaconda3/bin/python
from Bio import SeqIO
from Bio.Seq import Seq
import re
import Levenshtein  ###计算两个等长序列之间对应位置上不同字符的个数
from fuzzywuzzy import process
import os
##https://python.fasionchan.com/zh_CN/latest/practices/string-similarity.html
#def similarity(a, b):
#   return SequenceMatcher(None, a, b).ratio()
##https://blog.csdn.net/u010454729/article/details/46559845
dict1a={}
dict1b={}
dict2a={}
dict2b={}

# f1=open('test_r.fq','w')
# f2=open('test2_r.fq','w')
# f3=open('test3_r.fq','w')
f1=open('ATAC-110225_demultiplexed_R1.fq','w')
f2=open('ATAC-110225_demultiplexed_R2.fq','w')
fq1='ATAC-library-110225-ATAC-library-110225_FKDL202568736-1a-AK8856_1.clean.fq'
fq2='ATAC-library-110225-ATAC-library-110225_FKDL202568736-1a-AK8856_2.clean.fq'

print("############################## def loop_formacth1 ##############################")
def loop_formacth1(str1,longstr):
    dict_num={}
    for i in range(10,116):
        # print(longstr[i:(i+34)])
        # print(str1)
        sim1=Levenshtein.hamming(str1, longstr[i:(i+34)])
        dict_num[longstr[i:(i+34)]]=sim1
    dict_num=sorted(dict_num.items(), key=lambda x: x[1], reverse=False)
    # print(dict_num[0])
    return dict_num[0]
print("############################## def loop_formacth2 ##############################")
def loop_formacth2(str2,longstr):
    dict_num={}
    for i in range(10,117):
        sim2=Levenshtein.hamming(str2, longstr[i:(i+33)])
        dict_num[longstr[i:(i+33)]]=sim2
    dict_num=sorted(dict_num.items(), key=lambda x: x[1], reverse=False)
    # print(dict_num[0])
    return dict_num[0]

def loop_formacth3(str3,longstr):
    dict_num={}
    for i in range(0,(len(longstr)-8)):
        sim2=Levenshtein.hamming(str3, longstr[i:(i+8)])
        dict_num[longstr[i:(i+8)]]=sim2
    dict_num=sorted(dict_num.items(), key=lambda x: x[1], reverse=False)
    # print(dict_num[0])
    return dict_num[0]
def check_start(substr,allstr):
    pos=allstr.find(substr)
    if pos !=-1:
        return allstr[(pos+len(substr)):]
    else:
        return allstr
def check_startv2(substr,allstr,allstr1):
    pos=allstr.find(substr)
    if pos !=-1:
        return allstr1[(pos+len(substr)):]
    else:
        return allstr1

def DNA_reverse(sequence):
    sequence = sequence.upper()
    return sequence[::-1]

mismatch=5
# mismatch1=8
str1='CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
str2='AGCAGCCGTCGCAGTCTACACATATTCTCTGTC'
str2=DNA_reverse(str2) ###'CTGTCTCTTATACACATCTGACGCTGCCGACGA'
# str3='AGAGACAG'
# str4='GACAGAGA'
# str4=DNA_reverse(str4)
#print(str2_complement)
#print(str2)
bc='NNNNNNNNNNNNNNNN'
print("############################## seqio.parse  fq1 fq2##############################")
for seq_record in SeqIO.parse(fq1,"fastq"):
    # print(seq_record)
    tmp=seq_record.format("fastq").split('\n')
    dict1a[seq_record.id]=tmp[1]
    dict1b[seq_record.id]=tmp[3]
for seq_record in SeqIO.parse(fq2,"fastq"):
    tmp=seq_record.format("fastq").split('\n')
    dict2a[seq_record.id]=tmp[1]
    dict2b[seq_record.id]=tmp[3]
#aa='ssssssfdfhjjkkdsfdfgggf'
#vv='dfh'
#vv1='dsfd'
# pos=aa.find('fdfh')
# print(pos)
#m=re.search(vv+r'(.+)'+vv1[0:2],aa)
# if m:
    # print(m.groups()[0])
print( "############################## the last ##############################")
for i in dict1a.keys():
    (simstr1,sim1)=loop_formacth1(str1, dict1a[i])
    (simstr2,sim2)=loop_formacth2(str2, dict2a[i])
    if sim1<=mismatch and sim2<=mismatch:
        dict1a[i]=check_start('AGACAG',dict1a[i])
        dict2a[i]=check_start('AGACAG',dict2a[i])
        dict1b[i]=check_startv2('AGACAG',dict1a[i],dict1b[i])
        dict2b[i]=check_startv2('AGACAG',dict2a[i],dict2b[i])
        pos1=dict1a[i].find(simstr1)
        pos2=dict2a[i].find(simstr2)
        bc=dict2a[i][(pos2+33):(pos2+33+16)]
        if len(bc)==16:
            r1=dict1b[i][0:pos1]
            r2=dict2b[i][0:pos2]
            q1=dict1b[i][0:pos1]
            q2=dict2b[i][0:pos2]
            f1.write('@'+bc+':'+i+' 1'+'\n'+r1+'\n'+'+\n'+q1+'\n')
            f2.write('@'+bc+':'+i+' 2'+'\n'+r2+'\n'+'+\n'+q2+'\n')
        
print('##########################gzip###################################')
f1.close()
f2.close()
cmdstr1='gzip ATAC-110225_demultiplexed_R1.fq -f'
cmdstr2='gzip ATAC-110225_demultiplexed_R2.fq -f'
os.system(cmdstr1)
os.system(cmdstr2)
