#https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
import numpy as np
import pandas as pd
import scanpy
 
foloderlist = ['Hematopoietic','Proliferating_Hem','Osteoblast','Monocytes','SMC','Platelets','MSC']
rawoutputfolder = "./Dataset/Bone/GSE152285/Celltype/"

for folder in foloderlist:
    print(folder)
    outputfolder = rawoutputfolder+folder+'_only/'

# print("start loading")
# adata = scanpy.read_10x_mtx(outputfolder,
#                         var_names='gene_ids',
#                         cache=False)
# print(adata.shape)
# # adata= scanpy.read_text("C:/Users/whl19/Documents/Code/ctBuilderV2.0/dataset/singlecell/Bone/05regCARTIconHEAnm5_expression_mapped_superficial.txt",'\t',True )
# adata.write_csvs(outputfolder,sep='\t',skip_data=False)


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
                samplelist.append(line.strip().split(',')[0])

    print(len(samplelist))

    # load 
    print("write data")
    df = pd.read_csv(outputfolder+"/X.csv",delimiter=',', header=None, index_col=False) 
    print(df.shape)
    df = pd.DataFrame(df.values,  columns=genelist,index = samplelist )
    # df.reindex(samplelist)

    df.T.to_csv(outputfolder+"/GSE152285_"+folder+"_only_preprocessed_log.txt",sep='\t',index=True)