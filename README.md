# Methylome-of-clonal-ant_Obir-v5.4
# material: 
chr_Obir_v5.4_genomic.gff.txt (gff file with chromosome info)

0_GFF_to_temp_bed_13683.py

1_extract_single_isoform_from_gff_070821_v1.py

2_mk_whole_gene_gffbed_singlemRNA_070821_v3.py

3_mk_exon_intron_gffbed_singlemRNA_070821_v3.py

# Step0. mk gff.temp file for all genes (n = 13683)
0_GFF_to_temp_bed_13683.py

• input : chr_Obir_v5.4_genomic.gff

• output: chr_Obir_v5.4_genomic.gff.temp (includes all the isoforms, n = 13683)
  
# Step1. mk the list of single isofroms from gff file (n = 11868)
1_extract_single_isoform_from_gff_070821_v1.py

• input : chr_Obir_v5.4_genomic.gff

• output: Obir_maxlen_isoform_geneID_11868.txt
