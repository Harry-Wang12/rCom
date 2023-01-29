


import BioknowledgeOperate
import CommunicationChecking    
import ScoringOperate
import GeneExpressionOperate

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

    # CommunicationBoneCells = ScoringOperate.CalculateCoummunicationMatrixtofile(CRoutes,Bonecells,CLRweight = 0.75,CRRweight = 0.75, outputfile = Communication_bonesamplefile )

    # GeneExpressionOperate.writenamedsampletofile(CommunicationBoneCells,Communication_bonesamplefile)

    # print("11")
    # GeneExpressionDatasetfile = './Dataset/Bone/GSE152285/GSE152285_osteoblast_populationAveraged.txt'
    # RouteScoringFile = './Dataset/Bone/GSE152285/RouteScore_GSE152285_populationAveraged.txt'
    # CombinedCells = CommunicationChecking.AddingRouteinformationToSamples(GeneExpressionDatasetfile,RouteScoringFile)


    # LRBDfile = outputfile
    # outputfolder =  "./ResultAnalysis/GSE152285/"
    # LigandReceptorpairs = CommunicationChecking.IdentifiedLRCells(LRBDfile,CombinedCells,outputfolder,routeThreshold=0.1, geneThreshold=0.1)