

from cProfile import label
import pandas as pd


import scanpy as sc


studydataset = "Data4-GSE152285"
# studydataset = "Data3"

print("start!")
# inputfolder1 = './Dataset/Bone/GSE152285/Dll4(E12)/'
# inputfolder2 = './Dataset/Bone/GSE152285/Normal/'
outputfolder= './Dataset/'+studydataset+'/RAW_nonelog/'
outputfile = outputfolder+"/preprocessed.txt"


# normaladata = sc.read_10x_mtx('./Dataset/'+studydataset+'/Normal/',
#                         var_names='gene_symbols',
#                         cache=True)

# normaladata.var_names_make_unique()
# normaladata

# Dll4adata = sc.read_10x_mtx('./Dataset/'+studydataset+'/Dll4(E12)/',
#                         var_names='gene_symbols',
#                         cache=True)

# Dll4adata.var_names_make_unique()
# Dll4adata

adata = sc.read_10x_mtx('./Dataset/'+studydataset+'/raw/',
                        var_names='gene_symbols',
                        cache=True)

adata.var_names_make_unique()
# adata = normaladata
# adata = normaladata.concatenate(Dll4adata)

print(adata.shape)
# inputfolder = "./Dataset/Bone/Mouse bone study1/RAW/"
# outputfolder= './Dataset/Bone/Mouse bone study1/RAW_preocessed_new/'
# outputfile = outputfolder+"/preprocessed_log.txt"



# adata = sc.read_10x_mtx(inputfolder,
#                         var_names='gene_symbols',
#                         cache=False)

# adata.var_names_make_unique()

# adata = sc.read_csv('./Dataset/Bone/GSE104782/GSE104782_allcells_UMI_count.txt', delimiter='\t').transpose()



# labelfilepath = "./Dataset/"+studydataset+"/mc_compliant_metadata_example.tsv"
# samplecol=0
# labelcol =7

# # with open(labelfilepath,'r') as labelfile:
# #     for line in labelfile.readlines():
# #         samplename = line.strip().split('\t')[0]
# #         celltypedsamples.append(samplename)


# celltypelist={}
# with open(labelfilepath,"r") as celltypefile:
#     for line in celltypefile.readlines():
#         linedata= line.strip().split('\t')
        
#         celltypelist[linedata[samplecol]] = linedata[labelcol]


# orderedcelltype = []
# # print(adata.obs.keys())
# for x in adata.obs.index.values:
# #     print(x)
#     if x in celltypelist.keys():
#         orderedcelltype.append(celltypelist[x])
#     else:
#         orderedcelltype.append("unknown")
# adata.obs['celltype'] = orderedcelltype        
        
  
# celltypedsamples = []

# labelfilepath = "./Dataset/"+studydataset+"/meta_annotation.txt"
# sampenamecol = 0
# labelnamecol = 1

# with open(labelfilepath,'r') as labelfile:
#     for line in labelfile.readlines():
#         samplename = line.strip().split('\t')[0]
#         celltypedsamples.append(samplename)


# celltypelist={}
# with open('./Dataset/'+studydataset+'/meta_annotation.txt',"r") as celltypefile:
#     for line in celltypefile.readlines():
#         linedata= line.strip().split('\t')
#         celltypelist[linedata[0]] = linedata[-2]

# selectedcell = []
# orderedcelltype = []
# # print(adata.obs.keys())
# for x in adata.obs.index.values:
# #     print(x)
#     if x in celltypelist.keys():
#         orderedcelltype.append(celltypelist[x])
#         selectedcell.append(x)
#         print(len(selectedcell))
#     else:
#         orderedcelltype.append("unknown")
# adata.obs['celltype'] = orderedcelltype



print(adata.shape)

# adata = adata[selectedcell]


sc.pp.filter_cells(adata, min_genes=10)
sc.pp.filter_genes(adata, min_cells=100)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# sc.pp.normalize_total(adata, target_sum=10000)

# sc.pp.log1p(adata)

# sc.pp.scale(adata, max_value=10)


print(adata.shape)
adata.write_csvs(outputfolder, skip_data=False)


print("read genename ")
genelist = []
with open(outputfolder+"/var.csv",'r') as genenamefile:
    title = True
    for line in genenamefile.readlines():
        if title:
            title=False
        else:
            genelist.append(line.strip().split(',')[0])
print(len(genelist))


print("read samplename")
samplelist=[]
with open(outputfolder+"/obs.csv",'r') as samplenamefile:
    title = True
    for line in samplenamefile.readlines():
        if title:
            title=False
        else:
            samplename = line.strip().split(',')[0]

            samplelist.append(samplename)

print(len(samplelist))


# get labeled datset


# cellleidenlabel = GeneExpressionOperate.loadlabelfile(labelfilepath,sampenamecol,labelnamecol,"\t")


# load 
print("write data")
df = pd.read_csv(outputfolder+"/X.csv",delimiter=',', header=None, index_col=False) 
print(df.shape)
df = pd.DataFrame(df.values,  columns=genelist,index = samplelist )

# df.reindex(samplelist)
print(df.shape)

df.to_csv(outputfile,sep='\t',index=True)




# #########################################################################################
# sc.pp.filter_cells(adata, min_genes=5)
# sc.pp.filter_genes(adata, min_cells=10)



# adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
# sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)



# # adata = adata[adata.obs.n_genes_by_counts < 6000, :]
# # adata = adata[adata.obs.pct_counts_mt < 7.5, :]


# sc.pp.normalize_total(adata, target_sum=10000)
# # adata.write_csvs(outputfolder+"normalized/", skip_data=False)

# #Logarithmize the data:
# sc.pp.log1p(adata)

# sc.pp.scale(adata, max_value=10)

# print(adata.shape)
# adata.write_csvs(outputfolder, skip_data=False)


# print("read genename ")
# genelist = []
# with open(outputfolder+"/var.csv",'r') as genenamefile:
#     title = True
#     for line in genenamefile.readlines():
#         if title:
#             title=False
#         else:
#             genelist.append(line.strip().split(',')[0])
# print(len(genelist))


# print("read samplename")
# samplelist=[]
# with open(outputfolder+"/obs.csv",'r') as samplenamefile:
#     title = True
#     for line in samplenamefile.readlines():
#         if title:
#             title=False
#         else:
#             samplename = line.strip().split(',')[0]

#             samplelist.append(samplename)

# print(len(samplelist))


# # get labeled datset


# # cellleidenlabel = GeneExpressionOperate.loadlabelfile(labelfilepath,sampenamecol,labelnamecol,"\t")


# # load 
# print("write data")
# df = pd.read_csv(outputfolder+"/X.csv",delimiter=',', header=None, index_col=False) 
# print(df.shape)
# df = pd.DataFrame(df.values,  columns=genelist,index = samplelist )

# # df.reindex(samplelist)
# print(df.shape)

# df.to_csv(outputfile,sep='\t',index=True)



# #########################################################################################