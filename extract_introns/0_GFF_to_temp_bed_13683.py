import sys
#input : chr_Obir_v5.4_genomic.gff
#output: chr_Obir_v5.4_genomic.gff.temp > includes all the isoforms (n = 13683)


Input = ['chr_Obir_v5.4_genomic.gff']
#Input = ['repeat\\Obir.repeat.assembly.v5.4.fasta.out_tab_NoUnknown.gff']

print('Start:',Input[0])
outname = Input[0]+'.temp'
out = open(outname,'w')
with open(Input[0],'r') as fileopen:
    for i in fileopen.readlines():
        if i.startswith('Chr'):
            f = i.strip().split('\t')
            if f[2] == 'gene':
                print(f[0],f[3],f[4],f[8].split(';')[0].split('=gene-')[-1],f[5],f[6],file=out,sep='\t')
            else:
                pass
out.close()
            
print('end:',outname)
