

import GeneExpressionOperate


if __name__ =="__main__":
    dataset,genelist = GeneExpressionOperate.loadfilewithname( "./Dataset/Bone/GSE152285/GSE152285.csv",',')
    celllabel = GeneExpressionOperate.loadlabelfile("./Dataset/Bone/GSE152285/meta_annotation_revised.txt",0,-2)
    outputdir="./Dataset/Bone/GSE152285/"
    for celltypename,cellnamelist in celllabel.items():
        print(celltypename)
        subcells = GeneExpressionOperate.getsubgroupfromsamples(dataset,cellnamelist)
        outputfilepath = outputdir+celltypename+'.txt'
        GeneExpressionOperate.writenamedsampletofile(subcells,outputfilepath)
