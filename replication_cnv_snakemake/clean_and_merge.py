import pandas as pd
import numpy as np
import os
import argparse

## This dict refers to dataset PRJEB3235. Please change accordingly to your dataset
patients_dict = {1:'bc', 2:'bc',3:'bc',4:'bc',5:'bc',6:'bc',7:'bc',8:'sane',9:'sane',10:'sane',11:'sane',12:'sane',13:'sane',
                 14:'sane',15:'bc',16:'bc',17:'bc',18:'bc', 19:'bc', 20:'bc',21:'bc',22:'bc',23:'bc',24:'bc',25:'bc'}


parser = argparse.ArgumentParser(description='Merge results together.')

parser.add_argument('--prefix1', type=str, help= "Prefix of the files of the EXDE tool. Ex. \"ERR\", \"HI\"")
parser.add_argument('--prefix2', type=str, help= "Prefix of the files of the MOPS tool. Ex. \"ERR\", \"HI\"")
parser.add_argument('--paired_flag', type=int,help="0 if not paired end. 1 if it is.")
parser.add_argument('--result_folder',type=str, help="Folder where to save the results. Default: results")
parser.add_argument('--source_folder', type=str, help='Folder where to find the annotated files')
args = parser.parse_args()

prefix1 = args.prefix1
prefix2 = args.prefix2
paired_flag = args.paired_flag
result_folder = args.result_folder

if result_folder is None:
    result_folder = "results"

def concat_annotated_datasets(files, tool, qualityfilters = {}, lowcol = True):

    total = pd.DataFrame()
    for file in files:
        if 'annotated' not in file:
            continue
        data = pd.read_csv(file, sep="\t", low_memory=False)
        if 'HI' in file:
            data['record'] = file.split("/")[-1]

        total = pd.concat([total,data])

    total = total.loc[total['Annotation_mode'] == 'split']
    
    if tool == 'exde':
        total.rename(columns={'user#1':'start.p', 'user#2':'end.p', 'user#4':'nexons', 'user#9':'BF', 'user#10':'reads.expected', 
                              'user#11':'reads.observed', 'user#12':'reads.ratio'},inplace=True)
        if lowcol:
            total = total[['SV_chrom', 'SV_type', 'Gene_name', 'record']]

        
    if tool == 'cnmo':
        total['user#6'] = total['user#6'].str.replace('CN', '')
        total['user#6'] = total['user#6'].astype(int)
        total.loc[total['user#6']<2, 'SV_type'] = 'DEL'
        total.loc[total['user#6']>2, 'SV_type'] = 'DUP'
        total.rename(columns={'user#3':'record', 'user#4':'I/NI median', 'user#5':'I/NI mean', 'user#6':'CNV'},inplace=True)
        if lowcol:
            total = total[['SV_chrom', 'SV_type', 'Gene_name', 'record']]

    for key,value in qualityfilters.items():
        if key in total.columns:
            total = total.loc[total[key] >= value]

    def normalizestring(string):
        if '_' in string:
            string = "".join( string.split("_")[2:])
        string = string.split(".bam")[0]
        return string
    
    total['record'] = total['record'].apply(normalizestring)

    total.sort_values(by='record', inplace=True, ascending=True)

    i = 1
    d = 0
    hi = 19
    for patient in sorted( total.record.unique() ):
        if 'ERR' in patient:
            total.loc[total.record == patient,'patient'] = i
            if d == 1:
                i = i + 1
                d = 0
            else:
                d = d + 1
        
        if 'HI' in patient:
            total.loc[total.record == patient,'patient'] = hi
            hi = hi + 1


    total['patient'] = total['patient'].astype(int)

    total['disease'] = total['patient']
    total = total.replace({"disease": patients_dict})

    for col in total.columns:
        total[col] = total[col].astype(str)
        
    return total


def select_where_both_reads(df):
    test = df.groupby(['patient', 'SV_chrom', 'SV_type', 'Gene_name']).count().reset_index()
    test = test.loc[test['record'] > 1]
    test = test[['patient','SV_chrom','SV_type', 'Gene_name']]
    test = test.merge(df, on=['patient','SV_chrom','SV_type', 'Gene_name'], how = 'left')
    test = test.sort_values('BF', ascending=False) if 'BF' in df.columns else test.sort_values('I/NI median', ascending = False)
    test = test.drop_duplicates(subset=['patient','SV_chrom','SV_type','Gene_name'])

    return test


folder = os.path.join(args.source_folder)

exde_1 = [os.path.join(folder,f) for f in os.listdir(folder) if ( prefix1 in f )]
mops_1 = [os.path.join(folder,f)  for f in os.listdir(folder) if ( prefix2 in f )]


exde_1 = concat_annotated_datasets(exde_1, 'exde',lowcol=False)
mops_1 = concat_annotated_datasets(mops_1, 'cnmo',lowcol=False)


if paired_flag == 1:
    exde_1 = select_where_both_reads(exde_1)
    mops_1 = select_where_both_reads(mops_1)


filterlist_mops = ['patient', 'disease', 'AnnotSV_ID', 'Gene_name', 'SV_type', 'SV_chrom', 'Location', 'OMIM_ID', 'P_snvindel_nb' ,'P_snvindel_phen', 
                   'record', 'I/NI median', 'I/NI mean', 'CNV','B_gain_AFmax', 'B_loss_AFmax','B_ins_AFmax', 'B_inv_AFmax','AnnotSV_ranking_score', 'ACMG_class']
filterlist_exde = ['patient', 'disease', 'AnnotSV_ID', 'Gene_name', 'SV_type', 'SV_chrom', 'record', 'start.p', 'end.p', 'nexons', 'BF', 
                   'reads.expected', 'reads.observed', 'reads.ratio','B_gain_AFmax','B_loss_AFmax','B_ins_AFmax' ,'B_inv_AFmax','AnnotSV_ranking_score','ACMG_class']

filterlist_mops = [f for f in filterlist_mops if f in mops_1.columns]
filterlist_exde = [f for f in filterlist_exde if f in exde_1.columns]

if not os.path.exists(result_folder):
    os.makedirs(result_folder)

mops_1 = mops_1[filterlist_mops].reset_index(drop=True)
mops_1.to_excel(os.path.join(result_folder, f"mops_{prefix1}.xlsx") )
exde_1 = exde_1[filterlist_exde].reset_index(drop=True)
exde_1.to_excel(os.path.join(result_folder, f"exde_{prefix1}.xlsx"))

mergeon = ['SV_chrom','SV_type','Gene_name','record', 'patient', 'disease']
merge_1 = exde_1.merge(mops_1, on=mergeon, how='inner').drop_duplicates(subset=['SV_chrom','SV_type', 'Gene_name','patient'])


merge_1.rename(columns={'ACMG_class_x':'ACMG'}, inplace=True)
merge_1.drop(columns=['ACMG_class_y'], inplace=True)


merge_1['ACMG'] = merge_1['ACMG'].str.replace('full=', '')

merge_1.to_excel(os.path.join(result_folder, "total_1.xlsx"))