


import BioknowledgeOperate
# import CommunicationChecking    
import ScoringOperate
import GeneExpressionOperate
from multiprocessing import Process, Manager
import statistics
import os

def getnodeValue(nodetype,nodegenelist,sample):
    if nodetype == 'gene':
        if nodegenelist[0] in sample.keys():
            return sample[nodegenelist[0]]
        else:
            return 0

    elif nodetype =='And_bundle':
        effectivenumber = 0
        realvalue = 0
        subweight = 0
        for genename in nodegenelist:
            if genename in sample.keys():
                effectivenumber+=1
                realvalue+=sample[genename]
                if sample[genename]>0:
                    subweight+=1                    
                elif sample[genename]<0:
                    subweight-=1
        if effectivenumber==0:
            return 0
        else:
            return (abs(subweight)/effectivenumber)*(realvalue/effectivenumber)
    
    elif nodetype =='OR_bundle':
        effectivenumber = 0
        genevaluelist = []
        for genename in nodegenelist:
            if genename in sample.keys():
                effectivenumber+=1
                genevaluelist.append(sample[genename])
           
        if effectivenumber==0:
            return 0
        else:
            if max(genevaluelist)>0:
                return max(genevaluelist)               
            else:
                return statistics.mean(genevaluelist)  



def CLRcalculation(sample,CLR,weightratio = 0.75):
    # CLR = {        
    #     'Up':[],
    #     'Ligand':[],
    #     'Expectation':[],
    #     'Types':[]
    # }
    # the value suppose get value from 0-1
    Ligandvalue = getnodeValue(CLR['Types'][-1],CLR['Ligand'],sample)
    Upvalue = 0
    for i in range(len(CLR['Up'])):
        Upvalue+= CLR['Expectation'][i] * getnodeValue(CLR['Types'][i],CLR['Up'][i],sample)

    return weightratio*Ligandvalue+Upvalue*(1-weightratio)

def CRRcalculation(sample,CRR,weightratio = 0.75):
    # CRR = {
    #     'Receptor':[],        
    #     'Down':[],        
    #     'Expectation':[],
    #     'Types':[]        
    # }

    # the value suppose get value from 0-1

    Receptorvalue = getnodeValue(CRR['Types'][0],CRR['Receptor'],sample)
    Downvalue = 0
    for i in range(len(CRR['Down'])):
        Downvalue+= CRR['Expectation'][i+1]*getnodeValue(CRR['Types'][i+1],CRR['Down'][i],sample)
        
    return weightratio*Receptorvalue+Downvalue*(1-weightratio)



def getcellcommunicationVector(linelist,routeid,CommunicationRoute,samples):
    print("Start",routeid)
    CLRline = routeid+'-CLR'
    CRRline = routeid+'-CRR'
    for samplename,cellsample in samples.items():
        CLRscore = CLRcalculation(cellsample,CommunicationRoute['CLR'])
        CRRscore = CRRcalculation(cellsample,CommunicationRoute['CRR'])
        CLRline+='\t'+str(CLRscore)
        CRRline+='\t'+str(CRRscore)
    
    linelist.append(CLRline+'\n')
    linelist.append(CRRline+'\n')
    print("Done",routeid)
    





if __name__=='__main__':
    ChipSeqfilepath = "./DataBaseCombination/ChipSeqLibrary/gene_attribute_matrix.txt"

    ChipSeqDB= BioknowledgeOperate.loadChipSeqfile(ChipSeqfilepath)

    CellTalkfilepath = "./DataBaseCombination/CellTalkLRDB/human_lr_pair.txt"

    CellTalkDB = BioknowledgeOperate.loadCellTalkDB(CellTalkfilepath)

    CellChatfilepath = "./DataBaseCombination/CellChatLRDB/Human_complex_input_CellChatDB.csv"

    CellChatDB = BioknowledgeOperate.loadCellChatDB(CellChatfilepath)

    outputfile = "./DataBaseCombination/RouteStickedResult_V6.txt"

    rpacFile = "./DataBaseCombination/KEGG-Pathways/Kegg-Routes.txt"

    pathways_folder= "./DataBaseCombination/KEGG-Pathways/Raw/"

    bonesamplefile = "./Dataset/Bone/GSE152285/Raw_Normal_Processed/GSE152285_raw_Normal_preprocessed_log.txt"

    Communication_bonesamplefile = "./Dataset/Bone/GSE152285/Raw_Normal_Processed/Communication_GSE152285_raw_Normal_preprocessed_log.txt"

    iswrite = True

    socketPairs= BioknowledgeOperate.analysisSocketPairs(ChipSeqDB, CellTalkDB,CellChatDB)

    rpacRoutes,pathwaygraphmap = BioknowledgeOperate.splitrpacRoutes(rpacFile,pathways_folder,"")

    BioknowledgeOperate.generateCellCommunicationDB_V1(rpacRoutes,pathwaygraphmap,socketPairs,outputfile,iswrite)

    # CRoutes = ScoringOperate.loadCommunicationRoutes(outputfile,pathwaygraphmap)

    # Bonecells,genelist = GeneExpressionOperate.loadfilewithname(bonesamplefile,'\t')


    # writeallline = ScoringOperate.CalculateCoummunicationMatrixtolist(CRoutes,Bonecells,CLRweight = 0.75,CRRweight = 0.75 )
    
    # with open(Communication_bonesamplefile,'w') as Communication_bonesample: 
    #     for line in writeallline:
    #         Communication_bonesample.write(line)  

    # os.system("shutdown -s -t  1")

    # with Manager() as manager:
    #     newsamples=manager.dict()
    #     l = []
    #     for samplename,sample in Bonecells.items(): 
    #         index+=1
    #         # change to multiple process    
    #         p = Process(target=getcellcommunicationVector, args=(newsamples,samplename,sample, CRoutes,0.75 ,0.75,index))
    #         p.start()
    #     #     l.append(p)
        #     # p.join() 
        # for p in l:
        #     p.join()

    # GeneExpressionOperate.writenamedsampletofile(newsamples,Communication_bonesamplefile)