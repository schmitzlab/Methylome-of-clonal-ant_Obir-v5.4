# 1_cg_ratio_ant_gff_introns_062921_v3.2.py
description: calculate mCG ratio for each genomic contents such as;
'Genome','Genes','CDS','Exons','lnc_RNA','rRNA','snRNA','tRNA','TE','Introns',"5' UTR","3' UTR"
this code only considers the cytosines that has larger than 4 mapped reads (cut off >= 5)

usage: python3 1_cg_ratio_ant_gff_introns_062921_v3.2.py <'directory name'>

> make a directory and then put all of the <allc.tsv> files which need to be analyzed into that directory

input: multiple <allc.tsv> files in the specific directory

output: output1_Obir_genesTE_CG_ratio_co5_intron_070121.txt

# 2_methyl_Metaplot_v2_for_Ant_070621.py

usage: python3 2_methyl_Metaplot_v2_for_Ant_070621.py <gff.bed file> <size.genome file> <tsv.cw.bed file> <OUTPREFIX>

input: <> <ObirSize_chr.genome>
