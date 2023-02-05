import os, sys, glob, time
import numpy as np

#usage: python3 /home/hj43453/01.py_scripts/.....py
# 070821 only extract single isoform of mRNA 

'''
gene   =================  #confirmed that gene and mRNA are equal.
exon   ===   ===    ====
CDS      =   ===    ==     UTR = exon - CDS
intron    ===   ====       intron = gene - exon
'''
#######################################

def readgff(name):

    print(time.ctime())
    exonDic = {}  # {'id1':[exon1,exon2,...,], 'id2':[]..}   


    nSwitch = 0
    infile = open(name,'r')
    for i in infile.readlines():
        if i.startswith('Chr'):         
            spl = i.rstrip().split('\t')
            sChr,sGnomon,sClass,nS,nE,strand = spl[0],spl[1],spl[2],int(spl[3]),int(spl[4]),spl[6]   #Chr1 Gnomon	gene	1070	7288   +        
            #Gnomon: gene CDS exon lnc_RNA

            if sGnomon == 'Gnomon':
                if nSwitch == 1:
                    if sClass == 'mRNA':
                        nSwitch = 1


                    elif sClass == 'transcript' or sClass == 'lnc_RNA':
                        nSwitch = 0

                    elif sClass == 'exon':
                        
                        if sID.split('@')[0] in spl[8]:
                            sExonPos = spl[0]+'@'+spl[3]+'@'+spl[4] #Chr1@1070@1467
                            exonDic.setdefault(sID,[]).append(sExonPos)
                        else: pass
                            
                    elif sClass == 'CDS':
                        nSwitch = 0
                        
                else: #0
                    
                    if sClass == 'gene':
                        nSwitch = 1
                        sID = spl[8].split(';')[0].split('=gene-')[-1] +'@'+strand #LOC105288151@+
                        #rnaID = spl[8].split(';')[0].split('=rna-')[-1]#XM_026971221.1


                      
    print('exonDic',len(exonDic.keys()),time.ctime())
    infile.close()
    return exonDic      
 




        

def count_dic(exonDic,intronDic):

    n,m = 0,0
    for sKey in exonDic.keys():
        n += len(exonDic[sKey])
    for sKey in intronDic.keys():
        m += len(intronDic[sKey])            
    print('exonNum,intronNum:',n,m)



        



    
def test_gff(name):

    rnaDic = {}  # {'id1':[exon1,exon2,...,], 'id2':[]..}   

    rnaset = set()
    n,m = 0,0
    nSwitch = 0
    infile = open(name,'r')
    for i in infile.readlines():
        if i.startswith('Chr'):         
            spl = i.rstrip().split('\t')
            sChr,sGnomon,sClass,nS,nE,strand = spl[0],spl[1],spl[2],int(spl[3]),int(spl[4]),spl[6]   #Chr1 Gnomon	gene	1070	7288   +        
            #Gnomon: gene CDS exon lnc_RNA

            if sGnomon == 'Gnomon':

                if sClass == 'mRNA': #24156
                    rnaID  = spl[8].split(';')[0].split('=rna-')[1]
                    geneID = spl[8].split(';')[1].split('=gene-')[1]
                    rnalen = nE-nS
                           
                    rnaDic.setdefault(geneID,[]).append(rnaID + '@' + str(rnalen)) # {gene1:[r1@100,r2@200,r3@150]}
         
    return rnaDic


def select_rna(rnaDic):

    maxRNAdic = {}
    for sID in rnaDic.keys():
        temp = []
        for i in rnaDic[sID]:
            spl = i.split('@')
            temp.append(int(spl[1]))
            #end for
        x = np.array(temp)
        nMax = np.argmax(x)
        maxRNA = rnaDic[sID][nMax]
        maxRNAdic[sID] = maxRNA


    print(len(maxRNAdic))
    
    return maxRNAdic

  



######## MAIN #########################################################


#exonDic   = readgff('chr_Obir_v5.4_genomic.gff')


rnaDic = test_gff('chr_Obir_v5.4_genomic.gff')
maxRNAdic = select_rna(rnaDic)



outfile = open('Obir_maxlen_isoform_geneID_11868.txt','w')
for geneID in maxRNAdic.keys():
    rnaID = maxRNAdic[geneID].split('@')[0]
    print(rnaID,geneID,sep='\t',file=outfile)

outfile.close()


# 'gene' = 13339
# mRNA ID = 24156
# unique mRNA = 11868 ***


'''
gene 13339
mRNA 11418
lnc_RNA 1428
transcript 493'''
