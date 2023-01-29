import GeneExpressionOperate



if __name__=='__main__':


    foloderlist = ['SMC','Platelets','Proliferating_Hem','Hematopoietic','Osteoblast','Monocytes']
    Communication_dataset = {}
    
    for folder in foloderlist:
        print(folder)
        Communication_bonesamplefile = "./Dataset/Bone/GSE152285/Celltype/"+folder+".txt"

        Dataset,routelist = GeneExpressionOperate.loadfilewithname(Communication_bonesamplefile)

        for cellname,commumicationroutes in Dataset.items():
            Communication_dataset[cellname]= commumicationroutes
    
    outputfile = "./Dataset/Bone/GSE152285/Celltype/raw_GSE152285_preprocessed_log_0.5.txt"

    GeneExpressionOperate.writenamedsampletofile(Communication_dataset,outputfile)




























