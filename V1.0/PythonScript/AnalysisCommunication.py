




# 
from re import L


CommunicationfilePath = "./ResultAnalysis/GSE152285_04/Combination_67431.txt"

CellAnnotationPath = "./Dataset/Bone/GSE152285/meta_annotation_revised.txt"

outputfilePath = "./ResultAnalysis/GSE152285_04/Combination_67431_CellType.txt"

celltype={}

with open(CellAnnotationPath,'r') as annotation:
    title = True
    for line in annotation.readlines():
        if title:
            title = False
        else:
            linedata = line.strip().split('\t')
            sample = linedata[0]
            celllabel = linedata[-2]
            celltype[sample] = celllabel




with open(outputfilePath,'w') as output:
    with open(CommunicationfilePath,'r') as Communcation:
        for line in Communcation.readlines():
            
            
            linesingle = line.strip()
            linedata = linesingle.split('\t')
            # CommuncationTYPE = linedata[0]
            CommuncationCell = linedata[1]
            if CommuncationCell in celltype.keys():
                CommuncationCelltype = celltype[CommuncationCell]
                output.write( linesingle +'\t'+CommuncationCelltype+'\n')


















