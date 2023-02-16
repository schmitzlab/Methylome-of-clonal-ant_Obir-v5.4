import sys
import pandas as pd
import pybedtools as pbt
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time


Input_files   = ["chr_Obir_v5.4_genomic.gff.bed","ObirSize_chr.genome","test_allc_10.tsv.cw.bed",'tout'] #sys.argv[1:5] 
updown_win    = 1000  #up and down stream length
window_count  = 20    #window count
#CG_ylim       = 0.4
#CHG_ylim      = 0.2
#CHH_ylim      = 0.1
#All_ylim      = 0.4


modelist = ['CG','CHG','CHH']
vCG, vCHG, vCHH = modelist[0],modelist[1],modelist[2]


def making_bin(GFF,Input_files,updown_win,window_count,T_F):
    gff_bed       =  pbt.BedTool.from_dataframe(GFF)  #gff.bed
    up_bed        =  pbt.BedTool.flank(gff_bed,g=[size.genome],l=[1000bp window],r=0,s=T_F) 
    w_up_bed      =  pbt.BedTool.window_maker(up_bed,b=up_bed,n=window_count,reverse=T_F,i='srcwinnum') 

    down_bed      =  pbt.BedTool.flank(gff_bed,g=Input_files[1],l=0,r=updown_win,s=T_F)  #downstream 1kb
    w_down_bed    =  pbt.BedTool.window_maker(down_bed,b=down_bed,n=window_count,reverse=T_F,i='srcwinnum') #down 1kb > 10 bins (bed info)

    w_body_bed    =  pbt.BedTool.window_maker(gff_bed,b=gff_bed,n=window_count,reverse=T_F,i='srcwinnum') #body > 10 bins

    return w_up_bed,w_body_bed,w_down_bed


def make_bin_level(inter,`,Gene_id_ls):
    Gene_context_bin_dic = {}
    for i in Gene_id_ls:       #i = geneID
        for j in range(1,window_count+1):  #range(1,11) : j = 1,2,3,4,.. ,10
            for k in [vCG,vCHG,vCHH]:      # k = cg,chg,chh
                Gene_context_bin_dic[i+'_'+str(j)+'_'+k] = []   # {geneID_1_cg:[ratio1,ratio2,ratio3..], ..}

    for i in inter.itertuples():

        try:
            Gene_context_bin_dic[i[4]+'_'+i[15]].append(i[14]) #geneID_4_CWA = [ratio]
        except KeyError: pass
        try:
            Gene_context_bin_dic[i[4]+'_'+i[11]].append(i[14]) #geneID_4_CHG = [ratio]
        except KeyError: pass

    Bin_dic = {}
    for i in Gene_context_bin_dic.keys():  # i = LOC105283215_4_CG,  value = [ratio1,ratio2,ratio3,..]
        if len(Gene_context_bin_dic[i]) == 0:
            Bin_dic[i] = np.NaN
        else:
            Bin_dic[i] = sum(Gene_context_bin_dic[i])/len(Gene_context_bin_dic[i]) #Bin_dic = {geneID_4_cg:mean_ratio, ..}

    temp_CG_ls  = []
    temp_CHG_ls = []
    temp_CHH_ls = []
    for i in Gene_id_ls:     #temp_CG_ls = [ [gene1_1, _2, _3,.. _10],[gene2],[gene3],[],[],[],,,[] ]
        temp_CG_ls.append([])
        temp_CHG_ls.append([])
        temp_CHH_ls.append([])
        for j in range(1,window_count+1):  # j = 1,2,3,..,10
            temp_CG_ls[-1].append(Bin_dic[i+'_'+str(j)+'_'+vCG])   #  Bin_dic[geneID_j_CG] = mean_ratio 
            temp_CHG_ls[-1].append(Bin_dic[i+'_'+str(j)+'_'+vCHG])
            temp_CHH_ls[-1].append(Bin_dic[i+'_'+str(j)+'_'+vCHH])
                                                                    #category = up or bd or dw
    CG_df  = pd.DataFrame(temp_CG_ls, columns=[(category+'_'+str(x)) for x in range(1,window_count+1)], index = Gene_id_ls) 
    CHG_df = pd.DataFrame(temp_CHG_ls,columns=[(category+'_'+str(x)) for x in range(1,window_count+1)], index = Gene_id_ls)
    CHH_df = pd.DataFrame(temp_CHH_ls,columns=[(category+'_'+str(x)) for x in range(1,window_count+1)], index = Gene_id_ls)

    return CG_df,CHG_df,CHH_df
 
#################################################################### start
print(Input_files,time.ctime())

GFF        = pd.read_table(Input_files[0],header=None,names=['Chr','start','end','annotation','.','direction']) #gff.bed
GFF_plus   = GFF[GFF['direction'] == '+'] 
GFF_minus  = GFF[GFF['direction'] == '-'] 

Gene_id_ls = list(GFF.annotation) 


pw_up_bed,pw_body_bed,pw_down_bed =  making_bin(GFF_plus,Input_files,updown_win,window_count,False) #gff plus 
mw_up_bed,mw_body_bed,mw_down_bed =  making_bin(GFF_minus,Input_files,updown_win,window_count,True) #gff minus

print('Making bin',time.ctime())

w_up_bed   = pw_up_bed.cat(mw_up_bed,postmerge=False)     #plus minus 
w_body_bed = pw_body_bed.cat(mw_body_bed,postmerge=False) #plus minus 
w_down_bed = pw_down_bed.cat(mw_down_bed,postmerge=False) #plus minus 

print('Loding',time.ctime())

mC         = pd.read_table(Input_files[2],header=None)  
mC_bed     = pbt.BedTool.from_dataframe(mC)             

print('Locating mC',time.ctime())

inter_w_up_bed     = pbt.bedtool.BedTool.intersect(w_up_bed,mC_bed,wa=True,wb=True) 
inter_w_body_bed   = pbt.bedtool.BedTool.intersect(w_body_bed,mC_bed,wa=True,wb=True)
inter_w_down_bed   = pbt.bedtool.BedTool.intersect(w_down_bed,mC_bed,wa=True,wb=True)

table_w_up_inter   = pd.read_table(inter_w_up_bed.fn,header=None)  
table_w_body_inter = pd.read_table(inter_w_body_bed.fn,header=None)
table_w_down_inter = pd.read_table(inter_w_down_bed.fn,header=None)


Up_CG_df,Up_CHG_df,Up_CHH_df = make_bin_level(table_w_up_inter,  'up',Gene_id_ls) 
Bd_CG_df,Bd_CHG_df,Bd_CHH_df = make_bin_level(table_w_body_inter,'bd',Gene_id_ls) 
Dw_CG_df,Dw_CHG_df,Dw_CHH_df = make_bin_level(table_w_down_inter,'dw',Gene_id_ls)

CG_df  = pd.concat([Up_CG_df, Bd_CG_df ,Dw_CG_df ],axis=1)  
CHG_df = pd.concat([Up_CHG_df,Bd_CHG_df,Dw_CHG_df],axis=1)
CHH_df = pd.concat([Up_CHH_df,Bd_CHH_df,Dw_CHH_df],axis=1)



methyl_sum = pd.DataFrame([CG_df.mean(skipna=True),CHG_df.mean(skipna=True),CHH_df.mean(skipna=True),CG_df.count(),CHG_df.count(),CHH_df.count()],index=[vCG,vCHG,vCHH,'bincount','bincount','bincount']).T  #print including bincount
methyl_sum.to_csv(Input_files[3]+'_metaplot_summary.txt',sep='\t')




