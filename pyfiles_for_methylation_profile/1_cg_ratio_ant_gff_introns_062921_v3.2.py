import os, sys, glob, time

#usage: python3 /home/hj43453/01.py_scripts/.....py <'folder name'>
# input should be in the string form
# successfully make the set of gene,CDS,exon..
# cutoff >= 5
# add intron, UTR function
'''
gene   =================  #confirmed that gene and mRNA are equal.
exon   ===   ===    ====
CDS      =   ===    ==     UTR = exon - CDS
intron    ===   ====       intron = gene - exon
'''

sFolder = sys.argv[1] #'folder name'
#sFolder = os.getcwd() #print(sFolder) > this folder is the currently executed directory, not home dir

myList = glob.glob(sFolder+'/*chr.tsv')   
print(myList)

#######################################


def mk_gffdic():

    nGffDic = {}
    for i in range(14):
        sNum = str(i + 1)  #start from 1 to 14 (chr number)
        nGffDic['Chr'+sNum] = []
        for j in range(11): #number of classes: gene, cds, exon,...
            nGffDic['Chr'+sNum].append(set())

    return nGffDic

def mk_tempset(nS,nE):
    nTempSet = set()
    for num in range(nS,nE+1):
        nTempSet.add(num)
    return nTempSet

def readgff(name):

    nGffDic = mk_gffdic()  #make the dic below
    '''nGffDic = {'Chr1':[set0(gene),set1(CDS),set2(exon),set3(lnc_RNA),set4(rRNA),set5(snRNA),set6(tRNA),set7(TE),set8(intron),set9(5UTR),set10(3UTR)],'Chr2':[set(),set()..,'Chr14':[set(),set()]}'''
    exonDic = {}  # {'id1':[exon1,exon2,...,], 'id2':[]..}   
    #nClassDic = {'gene':[],'exon':[],'CDS':[],'lnc_RNA':[],'repeat':[]}
    infile = open(name,'r')
    for i in infile.readlines():
        if i.startswith('Chr'):        
            spl = i.rstrip().split('\t')
            sChr,sGnomon,sClass,nS,nE = spl[0],spl[1],spl[2],int(spl[3]),int(spl[4])   #Chr1 Gnomon	gene	1070	7288
            nSet = mk_tempset(nS,nE) 
            '''Gnomon: gene CDS exon lnc_RNA
               cmsearch: snRNA, rRNA, 
               tRNAscan-SE: tRNA
               RepeatMasker: TE'''

            if sGnomon == 'Gnomon': 
                if sClass == 'mRNA':
                    nGffDic[sChr][0].update(nSet)
                    sID = spl[8].split(';')[1].split('=gene-')[-1] #LOC105288151
                    rnaID = spl[8].split(';')[0].split('=rna-')[-1]#XM_026971221.1

                elif sClass == 'CDS':
                    nGffDic[sChr][1].update(nSet)

                elif sClass == 'exon':
                    
                    if rnaID in spl[8]:
                        nGffDic[sChr][2].update(nSet)
                        sExonPos = spl[0]+'@'+spl[3]+'@'+spl[4] #Chr1@1070@1467
                        exonDic.setdefault(sID,[]).append(sExonPos)
                    else: pass
                        #print('Error. not mRNA exon:',spl[8].split(';')[0])
                        

                elif sClass == 'lnc_RNA':
                    nGffDic[sChr][3].update(nSet)

            elif sGnomon == 'RepeatMasker':
                nGffDic[sChr][7].update(nSet)

            else:
                if sClass == 'rRNA': 
                    nGffDic[sChr][4].update(nSet)

                elif sClass == 'snRNA':
                    nGffDic[sChr][5].update(nSet)

                elif sClass == 'tRNA':
                    nGffDic[sChr][6].update(nSet)


    print('exonDic',len(exonDic.keys()))

    infile.close()
    print('.....complete reading gff file')
    return nGffDic,exonDic      
 

def mk_intron(exonDic):   #exonDic = {'id1':[exon1,exon2,...,], 'id2':[]..}
                                      #Chr1@1070@1467 #Chr1@1470@1568 #Chr1@1760@1856
    intronDic = {}                    #     i1= 1468,1469  i2 = 1569,1759
    for sID in exonDic.keys():
        exonList = exonDic[sID]
        nExonNum = len(exonList)  # exon list
        for num in range(nExonNum - 1): #num = 0,1,2,3..
            sChr = exonList[num].split('@')[0]
            n1exonEnd   = int(exonList[num].split('@')[2])   #1467
            n2exonStart = int(exonList[num+1].split('@')[1]) #1470
            
            sIntron = sChr + '@' + str(n1exonEnd+1) +'@'+ str(n2exonStart-1) #Chr1@1468@1469
            intronDic.setdefault(sID,[]).append(sIntron)

    print('intronDic',len(intronDic.keys()))   

    return intronDic


def mk_UTRset(nUTRset):

    '''nUTRset = (1,2,3,4, 25,26,27)'''
    


    return n5utrset,n3utrset


def mk_intronset(nGffDic,intronDic):
                                         
    '''nGffDic = {'Chr1':[set0(gene),set1(CDS),set2(exon),set3(lnc_RNA),set4(rRNA),set5(snRNA),set6(tRNA),set7(TE),set8(intron),set9(5UTR),set10(3UTR)],'Chr2':['''

    for geneID in intronDic.keys():
        dataList = intronDic[geneID]  # Chr1@1000@1300
        for data in dataList:
            f = data.split('@')
            sChr,nS,nE = f[0],int(f[1]),int(f[2])
            nSet = mk_tempset(nS,nE)
            nGffDic[sChr][8].update(nSet)

    for sChr in nGffDic.keys():  
        nUTRset = nGffDic[sChr][2] - nGffDic[sChr][1]                #UTR(set) = exon(set2) - cds(set1)
        nGffDic[sChr][9].update(nUTRset)
        '''n5utr,n3utr = mk_UTRset(nUTRset) ###
        #nGffDic[sChr][9].update(n5utr)
        #nGffDic[sChr][10].update(n3utr)'''

    return nGffDic


def countDNMT(name):

    print(name)
    infile = open(name,'r')
    '''nClassDic = {'gene':[],'CDS':[],'exon':[],'lnc_RNA':[],'rRNA':[],'snRNA':[],'tRNA':[],'TE':[]}''' 
    fClassList = [[],[],[],[],[],[],[],[],[],[],[]]  #number of lists = number of classes
    fGenomeList = []

    n,line = 0,0
    for i in infile.readlines():
        spl = i.rstrip().split('\t')
        sChr,context = spl[0],spl[3]
        nPosit = int(spl[1]) #cytosine residue
        methyl,total  = int(spl[4]),int(spl[5]) 
        
        n += 1
        if total >= 5:
            #line += 1
            if context[:2] == 'CG':
                fRatio = methyl/total
                fGenomeList.append(fRatio)

                for num in range(11): #class num.
                    try: 
                        nSet = nGffDic[sChr][num] #for each chromosome: [set0(gene),set1(CDS),set2(exon),set3(lnc_RNA),set4(rRNA),set5(snRNA),set6(tRNA),set7(TE)]
                    
                        if nPosit in nSet:
                            fClassList[num].append(fRatio)
                        else: pass
                    except KeyError: pass
                    #end of for
        #end of for

    fGenomeRatio = sum(fGenomeList)/len(fGenomeList) 

    fRatioList = []

    m = 0
    for myList in fClassList:
        m += 1
        if len(myList) == 0: 
            print(m,myList) 
            fRatioList.append(0)
        try:
            fRatioList.append(sum(myList)/len(myList))
        except ZeroDivisionError:
            pass

        #end of for
    
    infile.close()
    return fGenomeRatio,fRatioList
    print('process end....')



def calculate(file,fGenomeRatio,fRatioList):

    base = file.split('_L001_chr.')[0]
    sample = base.split('Kronauer_')[-1] #sample = D1_Cat_M2_S5
    
    print(sample,fGenomeRatio,fRatioList[0],fRatioList[1],fRatioList[2],fRatioList[3],fRatioList[4],fRatioList[5],fRatioList[6],fRatioList[7],fRatioList[8],fRatioList[9],fRatioList[10],sep='\t',file=out)



print(time.ctime())
out = open('Obir_genesTE_CG_ratio_co5_intron_070121.txt','w')
print('sample','Genome','Genes','CDS','Exons','lnc_RNA','rRNA','snRNA','tRNA','TE','Introns',"5' UTR","3' UTR",sep='\t',file=out)

nGffDic,exonDic = readgff('chr_Obir_v5.4_genomic.gff')  #gff_merge_Obir_v5.4_genes_TEnoUnkown.gff
print(time.ctime(),'....extracting introns')
intronDic = mk_intron(exonDic)
nGffDic = mk_intronset(nGffDic,intronDic)
print(time.ctime(),'gff over. start to read file')

for file in sorted(myList):
    fGenomeRatio,fRatioList = countDNMT(file)
    calculate(file,fGenomeRatio,fRatioList)
out.close()

print(time.ctime())



