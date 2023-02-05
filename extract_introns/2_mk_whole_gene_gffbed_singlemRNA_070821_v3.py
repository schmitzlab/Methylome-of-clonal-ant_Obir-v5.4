import os, sys, glob, time

#usage: python3 /home/hj43453/01.py_scripts/.....py
#this program makes gff.bed of the whole genes (maxlen single mRNA,n=11868)
#input: 'Obir_maxlen_isoform_geneID_11868.txt','chr_Obir_v5.4_genomic.gff.temp'
#output:chr_Obir_v5.4_genomic11868.gff.bed 

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
        rnaIDlist.append(i.rstrip().split('\t')[1]) #0 = rnaID, 1 = geneID
    return rnaIDlist 


def readgffbed(name,geneIDlist):

    outfile = open('chr_Obir_v5.4_genomic11868.gff.bed','w')


    n = 0
    infile = open(name,'r')
    for i in infile.readlines():
          
        spl = i.rstrip().split('\t')
        sID = spl[3]   #Chr1 Gnomon	gene	1070	7288   +        

        if sID in geneIDlist:
            print(i.rstrip(),file=outfile)

    outfile.close()


######## MAIN #########################################################


geneIDlist  = read_rnaID('Obir_maxlen_isoform_geneID_11868.txt')
readgffbed('chr_Obir_v5.4_genomic.gff.temp',geneIDlist)









