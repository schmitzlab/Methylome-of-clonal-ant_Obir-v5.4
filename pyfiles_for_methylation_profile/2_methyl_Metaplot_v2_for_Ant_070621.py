import sys,os
import pandas as pd
import pybedtools as pbt
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

#  information
#  bedtools 2.29 >
#  Input_files = [GENE ANNOTATION BED, GENOME BED, METHYLATION BED, OUTPREFIX]
#  'gffbed/exon1_Obirgene.gff.bed' < gff.bed (gene annotation bed) file should be located in the /gffbed/ directory.


Input_files   = sys.argv[1:5] 
updown_win    = 1000  #up and down stream length
window_count  = 20    #window bin
CG_ylim       = 0.4
CHG_ylim      = 0.2
CHH_ylim      = 0.1
All_ylim      = 0.4

args = sys.argv[5:]
option_dic={'-m':''}

if sys.argv[1:] == [] or sys.argv[1:] == ["-h"]:
    print(
"""==============================================
<DNA_methylation_metaplot Usage>
python3 *.py [Arapor11 gff.bed file] [genome file] [input tsv.cw.bed file] [out_file] -m <MODE>
Options
  -h           : help menu
  -m           : MODE option, 
                 -m CG (CG, CHG, CHH) < default 
                 -m CWA (CHH, CWA, nonCWA)
                 -m CWG (CHG, CWG, CCG)
==============================================""")    
    quit()

for i in range(len(args)):
    if args[i].startswith("-"):
        option_dic[args[i]]=args[i+1]
    else:
        pass#Input_files.append(args[i])

MODE = option_dic['-m'].upper()
#print(MODE)
if MODE == 'CG':
    modelist = ['CG','CHG','CHH']
elif MODE == 'CWG':  # CWG mode
    modelist = ['CHG','CWG','CCG']
    print('CWG mode')
elif MODE == 'CWA':  # CWA mode
    modelist = ['CHH','CWA','nonCWA']
    print('CWA mode')
else:
    print('MODE ERROR:',MODE)
    print('default mode is running: CG, CHG, CHH')
    modelist = ['CG','CHG','CHH']

vCG, vCHG, vCHH = modelist[0],modelist[1],modelist[2]
print(vCG, vCHG, vCHH)

def making_bin_updn(GFF,Input_files,updown_win,window_count,T_F):
    gff_bed       =  pbt.BedTool.from_dataframe(GFF)
    up_bed        =  pbt.BedTool.flank(gff_bed,g=Input_files[1],l=updown_win,r=0,s=T_F)
    w_up_bed      =  pbt.BedTool.window_maker(up_bed,b=up_bed,n=window_count,reverse=T_F,i='srcwinnum')

    down_bed      =  pbt.BedTool.flank(gff_bed,g=Input_files[1],l=0,r=updown_win,s=T_F)
    w_down_bed    =  pbt.BedTool.window_maker(down_bed,b=down_bed,n=window_count,reverse=T_F,i='srcwinnum')

    w_body_bed    =  pbt.BedTool.window_maker(gff_bed,b=gff_bed,n=window_count,reverse=T_F,i='srcwinnum')

    return w_up_bed,w_body_bed,w_down_bed
    #end

def making_bin_body(GFF,window_count,T_F):
    gff_bed       =  pbt.BedTool.from_dataframe(GFF)  #gff.bed
    w_body_bed    =  pbt.BedTool.window_maker(gff_bed,b=gff_bed,n=window_count,reverse=T_F,i='srcwinnum') #body > 10 bins
    return w_body_bed


def make_bin_level(inter,category,Gene_id_ls):
    Gene_context_bin_dic = {}
    for i in Gene_id_ls:
        for j in range(1,window_count+1):
            for k in [vCG,vCHG,vCHH]:
                Gene_context_bin_dic[i+'_'+str(j)+'_'+k] = []

    for i in inter.itertuples():
        #Index=3910, _1='Chr1', _2=24567417, _3=24567517, _4='LOC105281465_2', _5='Chr1', _6=24567425, _7=24567425, _8=11, _9='.', _10='-', _11='CG', _12=0, _13=11, _14=0.0, _15='.', _16='CGT')
        #Index=3910, _1='Chr1', _2=24567417, _3=24567517, _4='LOC105281465.4_2', _5='Chr1', 

        if '.' not in i[4]:  #_4='LOC105281465_2'
 
            try:
                Gene_context_bin_dic[i[4]+'_'+i[11]].append(i[14]) #'LOC105281465_4 _ CHG = [ratio]
            except KeyError: pass

        else:
            try:
                gene = i[4].split('_')[0].split('.')[0] + '_' + i[4].split('_')[1]  #_4='LOC105281465.4_2'
                Gene_context_bin_dic[gene+'_'+i[11]].append(i[14]) #geneID_4 _ CHG = [ratio]
            except KeyError: pass

    Bin_dic = {}
    for i in Gene_context_bin_dic.keys():
        if len(Gene_context_bin_dic[i]) == 0:
            Bin_dic[i] = np.NaN
        else:
            Bin_dic[i] = sum(Gene_context_bin_dic[i])/len(Gene_context_bin_dic[i])

    temp_CG_ls  = []
    temp_CHG_ls = []
    temp_CHH_ls = []
    for i in Gene_id_ls:
        temp_CG_ls.append([])
        temp_CHG_ls.append([])
        temp_CHH_ls.append([])
        for j in range(1,window_count+1):
            temp_CG_ls[-1].append(Bin_dic[i+'_'+str(j)+'_'+vCG])
            temp_CHG_ls[-1].append(Bin_dic[i+'_'+str(j)+'_'+vCHG])
            temp_CHH_ls[-1].append(Bin_dic[i+'_'+str(j)+'_'+vCHH])

    CG_df  = pd.DataFrame(temp_CG_ls, columns=[(category+'_'+str(x)) for x in range(1,window_count+1)], index = Gene_id_ls)
    CHG_df = pd.DataFrame(temp_CHG_ls,columns=[(category+'_'+str(x)) for x in range(1,window_count+1)], index = Gene_id_ls)
    CHH_df = pd.DataFrame(temp_CHH_ls,columns=[(category+'_'+str(x)) for x in range(1,window_count+1)], index = Gene_id_ls)

    return CG_df,CHG_df,CHH_df
    #end

def gff_read(gff_file):

    GFF        = pd.read_table(gff_file,header=None,names=['Chr','start','end','annotation','.','direction']) #gff.bed 
    gff_plus   = GFF[GFF['direction'] == '+'] 
    gff_minus  = GFF[GFF['direction'] == '-'] #

    Gene_id_ls = list(GFF.annotation) #
    return gff_plus,gff_minus,Gene_id_ls


def mk_windowbed(file):

    pgff_exon1,mgff_exon1,sExon1List = gff_read(file)
    exonNum = os.path.basename(file).split('_')[0]  #exon1
    
    pw_exon1 =  making_bin_body(pgff_exon1,window_count,False) #gff plus > 10 bin bed
    mw_exon1 =  making_bin_body(mgff_exon1,window_count,True)
    wbed_exon1 = pw_exon1.cat(mw_exon1,postmerge=False) #plus minus merge

    return exonNum, wbed_exon1, sExon1List


def mk_table(w_up_bed,mC_bed):

    inter_w_up_bed     = pbt.bedtool.BedTool.intersect(w_up_bed,mC_bed,wa=True,wb=True) #match tsv allC > gff.bin 
    table_w_up_inter   = pd.read_table(inter_w_up_bed.fn,header=None)  #
    return table_w_up_inter 

#################################################################### start


print(Input_files,time.ctime()) # Start

##0 read gff files(for upstream and dnstream bed)
gff_plus,gff_minus,Gene_id_ls = gff_read(Input_files[0])  #Input_files[0]: whole gene gff file
pw_up_bed,pw_genebody_bed,pw_down_bed =  making_bin_updn(gff_plus,Input_files,updown_win,window_count,False) #gff plus 
mw_up_bed,mw_genebody_bed,mw_down_bed =  making_bin_updn(gff_minus,Input_files,updown_win,window_count,True) #gff minus
w_up_bed = pw_up_bed.cat(mw_up_bed,postmerge=False)
w_down_bed = pw_down_bed.cat(mw_down_bed,postmerge=False)

##2 body - exon1,intron1, ...
exon1,wbed_exon1,sExonList1 = mk_windowbed('gffbed/exon1_Obirgene.gff.bed')  #exonNum = 'exon1',  #file = 'gffbed/exon1_Obirgene.gff.bed'
exon2,wbed_exon2,sExonList2 = mk_windowbed('gffbed/exon2_Obirgene.gff.bed')
exon3,wbed_exon3,sExonList3 = mk_windowbed('gffbed/exon3_Obirgene.gff.bed')
exon4,wbed_exon4,sExonList4 = mk_windowbed('gffbed/exon4_Obirgene.gff.bed')
exonN,wbed_exonN,sExonListN = mk_windowbed('gffbed/exonN_Obirgene.gff.bed')
intron1,wbed_intron1,sIntronList1 = mk_windowbed('gffbed/intron1_Obirgene.gff.bed')
intron2,wbed_intron2,sIntronList2 = mk_windowbed('gffbed/intron2_Obirgene.gff.bed')
intron3,wbed_intron3,sIntronList3 = mk_windowbed('gffbed/intron3_Obirgene.gff.bed')
intronN,wbed_intronN,sIntronListN = mk_windowbed('gffbed/intronN_Obirgene.gff.bed')
print('Making bin',time.ctime())

#print(w_up_bed)


##3 read and match mC tsv file
mC         = pd.read_table(Input_files[2],header=None)  # tsv.bed (allC )
mC_bed     = pbt.BedTool.from_dataframe(mC)             
print('Locating mC',time.ctime())

table_w_up_inter = mk_table(w_up_bed,mC_bed)
table_w_down_inter = mk_table(w_down_bed,mC_bed)

table_exon1 = mk_table(wbed_exon1,mC_bed)
table_exon2 = mk_table(wbed_exon2,mC_bed)
table_exon3 = mk_table(wbed_exon3,mC_bed)
table_exon4 = mk_table(wbed_exon4,mC_bed)
table_exonN = mk_table(wbed_exonN,mC_bed)
table_intron1 = mk_table(wbed_intron1,mC_bed)
table_intron2 = mk_table(wbed_intron2,mC_bed)
table_intron3 = mk_table(wbed_intron3,mC_bed)
table_intronN = mk_table(wbed_intronN,mC_bed)


##4 make bin level
print('cal methyl level',time.ctime())

Up_CG_df,Up_CHG_df,Up_CHH_df = make_bin_level(table_w_up_inter,  'up',Gene_id_ls) #input: matched and merged bed, up, gene list
Dw_CG_df,Dw_CHG_df,Dw_CHH_df = make_bin_level(table_w_down_inter,'dw',Gene_id_ls)

cg_e1,chg_e1,chh_e1 = make_bin_level(table_exon1,'e1',Gene_id_ls)
cg_e2,chg_e2,chh_e2 = make_bin_level(table_exon2,'e2',Gene_id_ls)
cg_e3,chg_e3,chh_e3 = make_bin_level(table_exon3,'e3',Gene_id_ls)
cg_e4,chg_e4,chh_e4 = make_bin_level(table_exon4,'e4',Gene_id_ls)
cg_eN,chg_eN,chh_eN = make_bin_level(table_exonN,'eN',Gene_id_ls)
cg_i1,chg_i1,chh_i1 = make_bin_level(table_intron1,'i1',Gene_id_ls)
cg_i2,chg_i2,chh_i2 = make_bin_level(table_intron2,'i2',Gene_id_ls)
cg_i3,chg_i3,chh_i3 = make_bin_level(table_intron3,'i3',Gene_id_ls)
cg_iN,chg_iN,chh_iN = make_bin_level(table_intronN,'iN',Gene_id_ls)



CG_df  = pd.concat([Up_CG_df,cg_e1,cg_i1,cg_e2,cg_i2,cg_e3,cg_i3,cg_e4,cg_iN,cg_eN,Dw_CG_df],axis=1)  #concat?  up,body,dn merge
CHG_df = pd.concat([Up_CG_df,cg_e1,cg_i1,cg_e2,cg_i2,cg_e3,cg_i3,cg_e4,cg_iN,cg_eN,Dw_CG_df],axis=1)
CHH_df = pd.concat([Up_CG_df,cg_e1,cg_i1,cg_e2,cg_i2,cg_e3,cg_i3,cg_e4,cg_iN,cg_eN,Dw_CG_df],axis=1)



print('making outfile',time.ctime())
methyl_sum = pd.DataFrame([CG_df.mean(skipna=True),CHG_df.mean(skipna=True),CHH_df.mean(skipna=True),CG_df.count(),CHG_df.count(),CHH_df.count()],index=[vCG,vCHG,vCHH,'bincount','bincount','bincount']).T  #print including bincount
methyl_sum.to_csv(Input_files[3]+'_metaplot_summary.txt',sep='\t')



print('end',time.ctime())




''' <CG_df>
              up_1  up_2  up_3  up_4  up_5  up_6  up_7  up_8  up_9  up_10  bd_1  bd_2  bd_3  bd_4  bd_5  bd_6  bd_7  bd_8  bd_9  bd_10  dw_1  dw_2  dw_3  dw_4  dw_5  dw_6  dw_7  dw_8  dw_9  dw_10
LOC105288151   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN
LOC105288134   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN
LOC113562166   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN
LOC105288138   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN
LOC105288140   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN
...            ...   ...   ...   ...   ...   ...   ...   ...   ...    ...   ...   ...   ...   ...   ...   ...   ...   ...   ...    ...   ...   ...   ...   ...   ...   ...   ...   ...   ...    ...
LOC105274780   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN
LOC105274781   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN
LOC105274810   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN
LOC105274784   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN
LOC105274785   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN    NaN
'''

''' <Up_CG_df>
LOC105287908_9_CHG nan
LOC105287908_9_CHH nan
LOC105287908_10_CG nan
LOC105287908_10_CHG nan
LOC105287908_10_CHH nan
LOC105287907_1_CG nan
LOC105287907_1_CHG nan
LOC105287907_1_CHH nan
LOC105287907_2_CG nan
LOC105287907_2_CHG nan
LOC105287907_2_CHH nan
'''

















