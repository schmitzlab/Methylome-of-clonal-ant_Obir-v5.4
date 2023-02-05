# Materials: 
chr_Obir_v5.4_genomic.gff.txt (gff file with chromosome info)

0_GFF_to_temp_bed_13683.py

1_extract_single_isoform_from_gff_070821_v1.py

2_mk_whole_gene_gffbed_singlemRNA_070821_v3.py

3_mk_exon_intron_gffbed_singlemRNA_070821_v3.py

# Step0. mk gff.temp file for all genes (n = 13683)
0_GFF_to_temp_bed_13683.py

• input : chr_Obir_v5.4_genomic.gff

• output: chr_Obir_v5.4_genomic.gff.temp (includes all the isoforms, n = 13683)
  
# Step1. mk the list of max length single isofroms (n = 11868) from the gff file 
1_extract_single_isoform_from_gff_070821_v1.py

• input : chr_Obir_v5.4_genomic.gff

• output: Obir_maxlen_isoform_geneID_11868.txt

# Step2. mk gff.bed file for all genes (n = 11868)
2_mk_whole_gene_gffbed_singlemRNA_070821_v3.py

• input: chr_Obir_v5.4_genomic.gff.temp (from step0), Obir_maxlen_isoform_geneID_11868.txt (from step1)

• output: chr_Obir_v5.4_genomic11868.gff.bed (this file will be used for metaplot)

# Step3. mk gff.bed files for exons and introns
3_mk_exon_intron_gffbed_singlemRNA_070821_v3.py

• input: chr_Obir_v5.4_genomic.gff, Obir_maxlen_isoform_geneID_11868.txt (from step1)

• output: gff.bed files for exon1,exon2,exon3,exon4,exonN,intron1,intron2,intron3,intronN
