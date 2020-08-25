#############猕猴的gtf转bed文件.py##########
print("#######################data input#####################################")
#f=open("100.bed","w")
f=open("Macaca_mulatta.Mmul_10.bed","w")
for line in open("Macaca_mulatta.Mmul_10.gtf"):
    # print("*******************spilt******************")
    lineL = line.split("\t")
    #print(lineL)
    # print('################染色体信息、判断基因类型##############')
    if len(lineL[0]) < 3:
        chrName = lineL[0]
        regionType = lineL[2]
        if (regionType == "gene"):
            # print('#############起始位置#################结束位置#')
            endgene = lineL[4]
            startgene = lineL[3]
            namegeneall = lineL[8]
            # print('################gtf最后一列genename、genesource....')
            # print(namegeneall)
            # print(type(namegeneall))
            y=namegeneall.split(";")
# gene_id "ENSMMUG00000000634"
# gene_version "4"
# gene_name "ZNF692"
# gene_source "ensembl"
# gene_biotype "protein_coding"
        # print(y[2])
        # print(y[2].split())
            l=y[2].split()
            c = "gene_name"
            if c == l[0]:
                # print('###################genename#######')
                genename=l[0]
                name=l[1]#################‘name’
                # print('#去掉name上的引号#')
                name=name.strip('"')
                #print(name)
                #print(chrName,startgene,endgene,genename,name)
                print("data write")
                f.write(chrName+'\t'+startgene+'\t'+endgene+'\t'+name+'\n')
f.close()
