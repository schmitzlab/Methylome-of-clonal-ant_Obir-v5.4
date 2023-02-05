import os, sys, glob, time

#usage: python3 /home/hj43453/01.py_scripts/.....py
# 070721 make all the gff.bed files for exons and introns
# input: 'Obir_maxlen_isoform_geneID_11868.txt','chr_Obir_v5.4_genomic.gff'
# output: gff.bed files for exon1,exon2,exon3,exon4,exonN,intron1,intron2,intron3,intronN

'''
gene   =================  #confirmed that gene and mRNA are equal.
exon   ===   ===    ====
CDS      =   ===    ==     UTR = exon - CDS
intron    ===   ====       intron = gene - exon
''' 
#######################################

def read_rnaID(name):

    rnaIDlist = []
    infile = open(name,'r')
    for i in infile.readlines():
        rnaIDlist.append(i.rstrip().split('\t')[0])
    return rnaIDlist 


def readgff(name,rnaIDlist):

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
                if sClass == 'exon':
                    if 'rna-X' in spl[8]:
                    
                        rnaID  = spl[8].split(';')[1].split('=rna-')[-1]
                        geneID = spl[8].split('gene=')[1].split(';')[0]
                        sID = geneID + '@' + strand                        #LOC105288151@+
                        if geneID.startswith('LOC'):pass
                        else: print('Error',geneID)

                        if rnaID in rnaIDlist:
                            sExonPos = spl[0]+'@'+spl[3]+'@'+spl[4] #Chr1@1070@1467
                            exonDic.setdefault(sID,[]).append(sExonPos)
                            #nSwitch += 1  #exon num = 82831
                        else: pass
                                        

    print('exonDic',len(exonDic.keys()),time.ctime()) #11868
    infile.close()
    return exonDic      
 

def mk_intronDic(exonDic):   #exonDic = {'id1@+':[exon1,exon2,...,], 'id2':[]..}
                                      #Chr1@1070@1467 #Chr1@1470@1568 #Chr1@1760@1856
    intronDic = {}                    #     i1= 1468,1469  i2 = 1569,1759
    for sID in exonDic.keys():
        strand = sID.split('@')[1] #+
        exonList = exonDic[sID]
        nExonNum = len(exonList)  # exon list
        for num in range(nExonNum - 1): #num = 0,1,2,3..

            if strand == '+':
                sChr = exonList[num].split('@')[0]
                n1exonEnd   = int(exonList[num].split('@')[2])   #1467
                n2exonStart = int(exonList[num+1].split('@')[1]) #1470
                sIntron = sChr + '@' + str(n1exonEnd+1) +'@'+ str(n2exonStart-1) #Chr1@1468@1469
                intronDic.setdefault(sID,[]).append(sIntron)
                
            elif strand == '-':   #  10@20   2@7
                sChr = exonList[num].split('@')[0]   
                n2exonEnd   = int(exonList[num+1].split('@')[2])   #7
                n1exonStart = int(exonList[num].split('@')[1]) #10               
                sIntron = sChr + '@' + str(n2exonEnd+1) +'@'+ str(n1exonStart-1) #Chr1@8@9
                intronDic.setdefault(sID,[]).append(sIntron)
            
            else: print('Error',sID,exonList)

            

    print('intronDic',len(intronDic.keys()))   

    return intronDic



    #exonDic = {'id1@+':[exon1,exon2,..], 'id2':[]..}
                      #Chr1@1070@1467 #Chr1@1470@1568 #Chr1@1760@1856
                      #     i1= 1468,1469  i2 = 1569,1759


    #intronDic = {'id1@+':[intron1,intron2, ..], 'id2':[]..}
                        #Chr1@1468@1469, Chr1@1520@1570, ...], 'id2':[]..}

                      

def mk_bedline(sKey,sInfo):  #Chr1@1070@1467 
    
    sID,strand = sKey.split('@')[0],sKey.split('@')[1] #id1, +
    spl = sInfo.split('@')
    sChr,sS,sE = spl[0],spl[1],spl[2]

    sLine = sChr +'\t'+ sS +'\t'+ sE +'\t'+ sID +'\t'+ '.' +'\t'+ strand


    return sLine   #Chr1    3069623	3073830	LOC113562264	.	-


def mk_otherline(sKey,myList,Num,outN):

    n = 0
    for sInfo in myList[Num:]:   #exonList[4:]
        n += 1
        sID,strand = sKey.split('@')[0],sKey.split('@')[1] #id1, +
        sID = sID + '.'+ str(n)
        
        spl = sInfo.split('@')
        sChr,sS,sE = spl[0],spl[1],spl[2]

        sLine = sChr +'\t'+ sS +'\t'+ sE +'\t'+ sID  +'\t'+ '.' +'\t'+ strand
        print(sLine,file = outN)



def mk_exonintron_bed(exonDic,intronDic,outList):

    for sKey in exonDic.keys():    #sKey = id1@+
        sExonList = exonDic[sKey] #[exon1,exon2,..]  exon1 = Chr1@1070@1467
        
        try:       
            sExon1 = mk_bedline(sKey,sExonList[0])  #extract 1st exon  #Chr1 3069623 3073830 LOC113562264 . -
            print(sExon1,file=outList[0]) #exon1)
            
            sExon2 = mk_bedline(sKey,sExonList[1])
            print(sExon2,file=outList[1]) #exon2)
            
            sExon3 = mk_bedline(sKey,sExonList[2])
            print(sExon3,file=outList[2]) #exon3)
            
            sExon4 = mk_bedline(sKey,sExonList[3])         
            print(sExon4,file=outList[3]) #exon4)
            
        except IndexError: pass
        
        try:
            mk_otherline(sKey,sExonList,4,outList[7])   #exonN
        except IndexError: pass

    for sKey in intronDic.keys():    #sKey = id1@+
        sIntronList = intronDic[sKey]
        
        try:
            sIntron1 = mk_bedline(sKey,sIntronList[0])
            print(sIntron1,file=outList[4]) #Intron1)
            
            sIntron2 = mk_bedline(sKey,sIntronList[1])
            print(sIntron2,file=outList[5]) #Intron2)
            
            sIntron3 = mk_bedline(sKey,sIntronList[2])
            print(sIntron3,file=outList[6]) #Intron3)
            
        except IndexError: pass
     
        try:
            mk_otherline(sKey,sIntronList,3,outList[8]) #intronN
        except IndexError: pass
        

def count_dic(exonDic,intronDic):

    n,m = 0,0
    for sKey in exonDic.keys():
        n += len(exonDic[sKey])
    for sKey in intronDic.keys():
        m += len(intronDic[sKey])            
    print('exonNum,intronNum:',n,m)

            
def outfile():
    
    exon1 = open('exon1_Obirgene.gff.bed','w')
    exon2 = open('exon2_Obirgene.gff.bed','w')
    exon3 = open('exon3_Obirgene.gff.bed','w')    
    exon4 = open('exon4_Obirgene.gff.bed','w')
    intron1 = open('intron1_Obirgene.gff.bed','w')
    intron2 = open('intron2_Obirgene.gff.bed','w')
    intron3 = open('intron3_Obirgene.gff.bed','w')
    exonN =  open('exonN_Obirgene.gff.bed','w')
    intronN = open('intronN_Obirgene.gff.bed','w')
    
    outList = [exon1,exon2,exon3,exon4,intron1,intron2,intron3,exonN,intronN]
            #    0       1    2     3      4        5       6     7      8
    return outList

def close_out(outList):
    for out in outList:
        out.close()


######## MAIN #########################################################

outList   = outfile() # [exon1,exon2,exon3,exon4,intron1,intron2,intron3,exonN,intronN]

rnaIDlist  = read_rnaID('Obir_maxlen_isoform_geneID_11868.txt')
exonDic   = readgff('chr_Obir_v5.4_genomic.gff',rnaIDlist)
intronDic = mk_intronDic(exonDic)
mk_exonintron_bed(exonDic,intronDic,outList)

count_dic(exonDic,intronDic)




close_out(outList)




