




import os




import GeneExpressionOperate



if __name__=='__main__':

    outputfolder= './Dataset/Bone/Mouse bone study1/Raw_processed/'
    celllogfilepath = outputfolder+"/preprocessed_log.txt"
     
    outputfilepath = outputfolder+"/radmonized_preprocessed_log.txt"


    GeneExpressionOperate.generaterandomdataset(celllogfilepath,outputfilepath)




    outputfolder= './Dataset/Bone/GSE152285/RAW/'
    celllogfilepath = outputfolder+"/GSE152285_preprocessed_log.txt"
     
    outputfilepath = outputfolder+"/radmonized_GSE152285_preprocessed_log.txt"


    GeneExpressionOperate.generaterandomdataset(celllogfilepath,outputfilepath)





    # os.system("shutdown -s -t  1")




