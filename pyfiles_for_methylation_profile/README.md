# 1_cg_ratio_ant_gff_introns_062921_v3.2.py
description: calculate mCG ratio for each genomic contents such as;
'Genome','Genes','CDS','Exons','lnc_RNA','rRNA','snRNA','tRNA','TE','Introns',"5' UTR","3' UTR"
this code only considers the cytosines that has larger than 4 mapped reads (cut off >= 5)

usage: python3 1_cg_ratio_ant_gff_introns_062921_v3.2.py <'directory name'>

> make a directory and then put all of the <sample_chr.tsv> files which need to be analyzed into that directory, and then run this code.

input: multiple <sample_chr.tsv> files in the specific directory

output: output1_Obir_genesTE_CG_ratio_co5_intron_070121.txt

# 2_methyl_Metaplot_v2_for_Ant_070621.py
description: calculate the average methylation ratio for each bin across the genes.

upstream of the gene start, each exon and intron, downstream of the gene end are splitted into 20 bins, and then average methyl ratio was calculated. 

before use: make 'gffbed' directory, and then put all the <exons.gff.bed> and <intons.gff.bed> into that directory.
the py file should be located right upstream of the 'gffbed' directory.

usage: python3 2_methyl_Metaplot_v2_for_Ant_070621.py <gff.bed file> <size.genome file> <tsv.cw.bed file> <OUTPREFIX>

usage example: python3 2_methyl_Metaplot_v2_for_Ant_070621.py <chr_Obir_v5.4_genomic11868.gff.bed> <ObirSize_chr.genome> <sample_chr.tsv.cw.bed> <D1cat_M3ei20bin>
  
*check 'extract_intron' folder to find the <gff.bed> example files

output: output2_Kro_D1cat_M3ei20bin_metaplot_summary.txt



# 3_DNA_methylation_Metaplot_cw_v2_github.py

description: metaplot for genes or repeat/TEs

same method with 2_methyl_Metaplot_v2_for_Ant_070621.py
  
  ##use <chr_Obir_v5.4_genomic11868.gff.bed> for gene analysis
  
  ##use <over100.Obir.repeat.NoUnknown.gff.bed> for repeat/TE analysis
