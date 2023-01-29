
import numpy as np
import pandas as pd
import scanpy as sc





adata1 = sc.read_csv('E:/Documents/codes/BOSC2022/Dataset/Good/Singlecell/GSE98969 AD_WT 2017/GSE98969_log2 AD1m_WT1m_AVG/ADvsControl/normalized_AD1mvsWT1m_out_all.csv',delimiter=',').transpose()

adata1.var_names_make_unique()



adata2 = sc.read_csv('E:/Documents/codes/BOSC2022/Dataset/Good/Singlecell/GSE98969 AD_WT 2017/GSE98969_log2 AD3m_WT3m_AVG/ADvsControl/normalized_AD3mWT3m_group_out_all.csv',delimiter=',').transpose()

adata2.var_names_make_unique()


adata3 = sc.read_csv('E:/Documents/codes/BOSC2022/Dataset/Good/Singlecell/GSE98969 AD_WT 2017/GSE98969_log2 AD6m_WT6m_AVG/ADvsControl/normalized_AD6mWT6m_group_out_all.csv',delimiter=',').transpose()

adata3.var_names_make_unique()




sc.pp.log1p(adata1)
sc.pp.log1p(adata2)
sc.pp.log1p(adata3)

sc.pp.scale(adata1, max_value=10)
adata1.write_csvs("E:/Documents/codes/BOSC2022/Dataset/Good/Singlecell/ArtificalMade/AD_1months_AD_WT_mouse/", skip_data=False)
sc.pp.scale(adata2, max_value=10)
adata2.write_csvs("E:/Documents/codes/BOSC2022/Dataset/Good/Singlecell/ArtificalMade/AD_3months_AD_WT_mouse/", skip_data=False)
sc.pp.scale(adata3, max_value=10)
adata3.write_csvs("E:/Documents/codes/BOSC2022/Dataset/Good/Singlecell/ArtificalMade/AD_6months_AD_WT_mouse/", skip_data=False)










