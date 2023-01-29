

from rPAC import rpac_route as route , rpac_score as score,rpac_pathway_graph as pathway_graph,rpac_util as util
from multiprocessing import Process, Manager
import BioknowledgeOperate
# import GeneExpressionOperate
# import statistics
# import psutil
import gc
import os
import pandas as pd

def SplitCommunicationRoute(linedata,G1,G2=False):
    ResultsHeader = ['id','C2-R2-ID','C2-R2-PW','C2-R2-TYPE','C2-R2-ids','C2-R2-genes','C2-R2-RRecetors-ids','C2-R2-RRecetors-genes','C2-R2-RReceptors-type','C2-R2-RLigands-ids','C2-R2-RLigands-genes','C2-R2-RLigands-type','SocketPairs_Receptor','SocketPairs_Ligand','SocektPairs_TF','C1-R1-TFs-ids','C1-R1-TFs-genes','C1-R1-ids','C1-R1-genes','C1-R1-PW','C1-R1-ID','C1-R1-TYPE','combineType']

    
    CLR = {
        
        'Up':[],
        'Ligand':[],
        'Expectation':[],
        'Types':[]
    }

    CRR = {
        'Receptor':[],        
        'Down':[],        
        'Expectation':[],
        'Types':[]
        
    }


    # get CLR
    if  linedata[ResultsHeader.index('combineType')] == '4':
        # get CLR
        CLR['Up'].append([linedata[ResultsHeader.index('SocektPairs_TF')]])
        CLR['Ligand'].append(linedata[ResultsHeader.index('SocketPairs_Ligand')])
        CLR['Expectation'].append(1)
        CLR['Expectation'].append(1)
        CLR['Types'].append('gene')
        CLR['Types'].append('gene')

        # get CRR
        CRR['Receptor'].append(linedata[ResultsHeader.index('SocketPairs_Receptor')])
        CRRdownnodes = linedata[ResultsHeader.index('C2-R2-ids')].split(',')
        CRRRids = linedata[ResultsHeader.index('C2-R2-RRecetors-ids')].split(',')
        CRRRexcatid = ""
        for id in CRRRids:
            if len(list(set(CRR['Receptor']).intersection(BioknowledgeOperate.getnodegenename(id,G1))))>0:
                CRRRexcatid = id
        addtodown = False
        for id in CRRdownnodes:
            if addtodown:
                Expectation = score.get_node_expectation(G1, previousid, id, CRR['Expectation'][-1])
                CRR['Down'].append(BioknowledgeOperate.getnodegenename(id,G1))
                CRR['Types'].append(BioknowledgeOperate.getnodetype(id,G1))
                CRR['Expectation'].append(Expectation)
                previousid= id
            elif id == CRRRexcatid:
                addtodown=True
                CRR['Types'].append('gene')
                previousid = id
                Expectation = 1
                CRR['Expectation'].append(Expectation)

    if  linedata[ResultsHeader.index('combineType')] == '5':   
        # get CLR
        # CLR['Up'].append(linedata[ResultsHeader.index('SocektPairs_TF')])
        CLR['Ligand'].append(linedata[ResultsHeader.index('SocketPairs_Ligand')])
        
        CLR['Expectation'].append(1)
        
        CLR['Types'].append('gene')

        # get CRR
        CRR['Receptor'].append(linedata[ResultsHeader.index('SocketPairs_Receptor')])
        CRRdownnodes = linedata[ResultsHeader.index('C2-R2-ids')].split(',')
        CRRRids = linedata[ResultsHeader.index('C2-R2-RRecetors-ids')].split(',')
        CRRRexcatid = ""
        for id in CRRRids:
            if len(list(set(CRR['Receptor']).intersection(BioknowledgeOperate.getnodegenename(id,G1))))>0:
            
                CRRRexcatid = id
        addtodown = False
        for id in CRRdownnodes:
            if addtodown:
                Expectation = score.get_node_expectation(G1, previousid, id, CRR['Expectation'][-1])
                CRR['Down'].append(BioknowledgeOperate.getnodegenename(id,G1))
                CRR['Types'].append(BioknowledgeOperate.getnodetype(id,G1))
                CRR['Expectation'].append(Expectation)
                previousid = id
            elif id == CRRRexcatid:
                addtodown=True
                previousid = id
                CRR['Types'].append('gene')
                Expectation = 1
                CRR['Expectation'].append(Expectation) 
        
    if  linedata[ResultsHeader.index('combineType')] == '2':   
        # get CLR
        # CLR['Up'].append(linedata[ResultsHeader.index('SocektPairs_TF')])
        CLR['Ligand']=linedata[ResultsHeader.index('SocketPairs_Ligand')].split(',')
        CLRuppernodes = linedata[ResultsHeader.index('C1-R1-ids')].split(',')
        
        CLRLexcatid = ""
        for id in CLRuppernodes:
            if len(list(set(CLR['Ligand']).intersection(BioknowledgeOperate.getnodegenename(id,G2))))>0:
            
                CLRLexcatid = id

        addtoup = False
        for id in CLRuppernodes[::-1]:
            if addtoup:
                Expectation = score.get_node_expectation(G2,  id,postid, CLR['Expectation'][0])
                CLR['UP'].insert(0,BioknowledgeOperate.getnodegenename(id,G2))
                CLR['Types'].insert(0,BioknowledgeOperate.getnodetype(id,G2))
                CLR['Expectation'].insert(0,Expectation)
                postid = id
            elif id == CLRLexcatid:
                addtodown=True
                postid = id
                Expectation = 1
                CLR['Types'].append(BioknowledgeOperate.getnodetype(id,G2))
                CLR['Expectation'].append(Expectation) 

        # get CRR
        CRR['Receptor']=linedata[ResultsHeader.index('SocketPairs_Receptor')].split(',')
        CRRdownnodes = linedata[ResultsHeader.index('C2-R2-ids')].split(',')
        CRRRids = linedata[ResultsHeader.index('C2-R2-RRecetors-ids')].split(',')
        CRRRexcatid = ""
        for id in CRRRids:
            if len(list(set(CRR['Receptor']).intersection(BioknowledgeOperate.getnodegenename(id,G1))))>0:
            # if CRR['Receptor'] in BioknowledgeOperate.getnodegenename(id,G1):
                CRRRexcatid = id
        addtodown = False
        for id in CRRdownnodes:
            if addtodown:
                Expectation = score.get_node_expectation(G1, previousid, id, CRR['Expectation'][-1])
                CRR['Down'].append(BioknowledgeOperate.getnodegenename(id,G1))
                CRR['Types'].append(BioknowledgeOperate.getnodetype(id,G1))
                CRR['Expectation'].append(Expectation)
                previousid = id
            elif id == CRRRexcatid:
                addtodown=True
                previousid = id
                Expectation = 1
                CRR['Types'].append('gene')
                CRR['Expectation'].append(Expectation)     
    return CLR,CRR

def loadCommunicationRoutes(CommunicationRoutefilepath,pathwaygraphmap):

    ResultsHeader = ['id','C2-R2-ID','C2-R2-PW','C2-R2-TYPE','C2-R2-ids','C2-R2-genes','C2-R2-RRecetors-ids','C2-R2-RRecetors-genes','C2-R2-RReceptors-type','C2-R2-RLigands-ids','C2-R2-RLigands-genes','C2-R2-RLigands-type','SocketPairs_Receptor','SocketPairs_Ligand','SocektPairs_TF','C1-R1-TFs-ids','C1-R1-TFs-genes','C1-R1-ids','C1-R1-genes','C1-R1-PW','C1-R1-ID','C1-R1-TYPE','combineType']

    CRoutes = {}
    RelatedGenes=[]



    with open(CommunicationRoutefilepath,'r') as Communicationfile:
        for line in Communicationfile.readlines()[1:]:
            linedata = line.strip().split('\t')
            Routeid = 'COMMUNICATIONROUEID_'+linedata[ResultsHeader.index('id')].upper()
            
            # if(linedata[ResultsHeader.index('id')]=='7789'):
            #     print('!')
            pw1_name =linedata[ResultsHeader.index('C2-R2-PW')]
            pw2_name = linedata[ResultsHeader.index('C1-R1-PW')]
            if not pw1_name=='NA':                
                G1 = pathwaygraphmap[pw1_name]
            if not pw2_name=='NA':
                G2 = pathwaygraphmap[pw2_name]
            else:
                G2 =False
            CLR,CRR = SplitCommunicationRoute(linedata,G1,G2)
            CRoutes[Routeid]={
                'CLR':CLR,
                'CRR':CRR
            }
            for genelist in CLR['Up']:
                for gene in genelist:
                    if gene not in RelatedGenes:
                        RelatedGenes.append(gene)
            for gene in CLR['Ligand']:
                if gene not in RelatedGenes:
                        RelatedGenes.append(gene)
            
            for genelist in CRR['Down']:
                for gene in genelist:
                    if gene not in RelatedGenes:
                        RelatedGenes.append(gene)
            for gene in CRR['Receptor']:
                if gene not in RelatedGenes:
                        RelatedGenes.append(gene)
        

        print("Find "+ linedata[ResultsHeader.index('id')]+" Communication Routes")
    return CRoutes,RelatedGenes

def getRouteline_multiple(writeingline,sample,genelist, CommunicationRoutes,CLRweight ,CRRweight,Cindex,totalnumber):
    
    print(str(Cindex)+'/'+str(totalnumber))    
    samplelist = []
    
    # genelist = genedict.keys()
    for routename,CommunicationRoute in CommunicationRoutes.items():
        # print(routename)
        # if routename == "CommunicationRoueid_6612":
        #     print("xx")
        # CLR
        # 
        # 
    #    getnodeValue_mulitple(nodetype,nodegenelist,sample,expectation,genelist)
        Ligandvalue = getnodeValue_mulitple(CommunicationRoute['CLR']['Types'][-1],CommunicationRoute['CLR']['Ligand'],sample,1,genelist)
        if Ligandvalue >0:
            Upvalue = 0
            if len(CommunicationRoute['CLR']['Up'])>0:
                for i in range(len(CommunicationRoute['CLR']['Up'])):
                    Upvalue +=  getnodeValue_mulitple(CommunicationRoute['CLR']['Types'][i],CommunicationRoute['CLR']['Up'][i],sample,CommunicationRoute['CLR']['Expectation'][i],genelist)
                
                samplelist.append (str(round(CLRweight*Ligandvalue+(Upvalue/(len(CommunicationRoute['CLR']['Types'])-1))*(1-CLRweight),3)))
            else:
                samplelist.append (str(Ligandvalue))
        else:
            samplelist.append (str(0))


        # CLRscore =CLRcalculation(cellsample,CommunicationRoute['CLR'],CLRweight)
        
        # CRR

        Receptorvalue = getnodeValue_mulitple(CommunicationRoute['CRR']['Types'][0],CommunicationRoute['CRR']['Receptor'],sample,1,genelist)
        if Receptorvalue >0:
            Downvalue = 0
            if len(CommunicationRoute['CRR']['Down'])>0:
                for i in range(len(CommunicationRoute['CRR']['Down'])):
                    Downvalue+= getnodeValue_mulitple(CommunicationRoute['CRR']['Types'][i+1],CommunicationRoute['CRR']['Down'][i],sample,CommunicationRoute['CRR']['Expectation'][i+1],genelist)
                samplelist.append(str(round(CRRweight*Receptorvalue+(Downvalue/(len(CommunicationRoute['CRR']['Types'])-1))*(1-CRRweight),3)))
            else:
                samplelist.append(str(Receptorvalue))
        else:
            samplelist.append (str(0))

            
        


        # CRRscore = CRRcalculation(cellsample,CommunicationRoute['CRR'],CRRweight)

    
    # resultlist.append(CRRline+'\n')
    print("Done "+str(Cindex)+'/'+str(totalnumber))
    writeingline.append(samplelist)
    # return samplelist
    # newsamples[samplename] = newsample
    # print("Done with "+ samplename)
    # return newsample

def getnodeValue_mulitple_v2(nodetype,nodegenelist,sample,expectation,genelist):

    if nodetype == 'gene':
        # if expectation>0:
        if nodegenelist[0] in genelist:
            Locindex = genelist.index(nodegenelist[0])
            
            return expectation*sample[Locindex]
            
        else:
                return 0

        # if expectation<0:
        #     if nodegenelist[0] in genelist:
        #         Locindex = genelist.index(nodegenelist[0])
        #         if sample[Locindex]<0:
        #             return abs(sample[Locindex])
        #         else:
        #             return 0
        #     else:
        #             return 0
        

    elif nodetype =='And_bundle':
        
        realvaluelist = [0]


        for genename in nodegenelist:
            if genename in genelist: 
                Locindex = genelist.index(genename)  
                
                realvaluelist.append(sample[Locindex])

        nodevalue = min(realvaluelist)
        return expectation*(nodevalue)



    elif nodetype =='OR_bundle':
        realvaluelist = [0]  
     
        for genename in nodegenelist:
            if genename in genelist:
                Locindex = genelist.index(genename) 
                realvaluelist.append(sample[Locindex])

        nodevalue = max(realvaluelist)

        return expectation *(nodevalue)

def getnodeValue_mulitple(nodetype,nodegenelist,sample,expectation,genelist):

    if nodetype == 'gene':
        if expectation>0:
            if nodegenelist[0] in genelist:
                Locindex = genelist.index(nodegenelist[0])
                if sample[Locindex]>0:
                    return sample[Locindex]
                else:
                    return 0
            else:
                    return 0

        if expectation<0:
            if nodegenelist[0] in genelist:
                Locindex = genelist.index(nodegenelist[0])
                if sample[Locindex]<0:
                    return sample[Locindex]
                else:
                    return 0
            else:
                    return 0
        

    elif nodetype =='And_bundle':
        
        realvaluelist = [0]
        getzero = False

        for genename in nodegenelist:
            if genename in genelist: 
                Locindex = genelist.index(genename)  
                if expectation > 0 and sample[Locindex]<0:
                    getzero=True
                realvaluelist.append(sample[Locindex])

        nodevalue = min(realvaluelist)

        if getzero:
            return 0
        else:
            return nodevalue



    elif nodetype =='OR_bundle':
        realvaluelist = [0]
        getzero = False
     
        for genename in nodegenelist:
            if genename in genelist:
                Locindex = genelist.index(genename) 

                if expectation < 0 and sample[Locindex]>0:
                    getzero=True
                realvaluelist.append(sample[Locindex])

        nodevalue = max(realvaluelist)

        if getzero:
            return 0
        else:
            return nodevalue



if __name__=="__main__":
    # load communication routes


    manager = Manager() 
    # ChipSeqfilepath = "./DataBaseCombination/ChipSeqLibrary/gene_attribute_matrix.txt"

    # ChipSeqDB= BioknowledgeOperate.loadChipSeqfile(ChipSeqfilepath)
    # print('RAM memory % used:', psutil.virtual_memory()[2])

    # # CellTalkfilepath = "./DataBaseCombination/CellTalkLRDB/human_lr_pair.txt"

    # # CellTalkDB = BioknowledgeOperate.loadCellTalkDB(CellTalkfilepath)
    # print('RAM memory % used:', psutil.virtual_memory()[2])

    # # CellChatfilepath = "./DataBaseCombination/CellChatLRDB/Human_complex_input_CellChatDB.csv"

    # # CellChatDB = BioknowledgeOperate.loadCellChatDB(CellChatfilepath)
    # print('RAM memory % used:', psutil.virtual_memory()[2])

    communicationfilepath = "./DataBaseCombination/RouteStickedResult_V6.txt"

    rpacFile = "./DataBaseCombination/KEGG-Pathways/Kegg-Routes.txt"

    pathways_folder= "./DataBaseCombination/KEGG-Pathways/Raw/"

    rpacRoutes,pathwaygraphmap = BioknowledgeOperate.splitrpacRoutes(rpacFile,pathways_folder,"")

    CRoutes,RouterelatedGene = loadCommunicationRoutes(communicationfilepath,pathwaygraphmap)

    # load file and create multiple thread

    # labelfilepath = "./Dataset/Bone/GSE152285/meta_annotation_revised.txt"
    # sampenamecol = 0
    # labelnamecol = 1

    # celltypedsamples = {}

    # with open(labelfilepath,'r') as labelfile:
    #     for line in labelfile.readlines():
    #         # samplename = line.strip().split('\t')[0]
    #         celltypedsamples[line.strip().split('\t')[sampenamecol]] = line.strip().split('\t')[labelnamecol]

   


    outputfolder= './Dataset/Bone/GSE152285/RAW_test6/'
    celllogfilepath = outputfolder+"/preprocessed.txt"
     
    Communication_bonesamplefile = './Dataset/Bone/GSE152285/RAW/Communication_preprocessed.txt'

    # newCroutes = {}

    # newCroutes["COMMUNICATIONROUEID_605"] = CRoutes["COMMUNICATIONROUEID_605"]
    # newCroutes["COMMUNICATIONROUEID_96376"] = CRoutes["COMMUNICATIONROUEID_96376"]

    # CRoutes = newCroutes

    # outputfolder= './Dataset/Bone/Mouse bone study1/RAW_unlog/'
    # # celllogfilepath = outputfolder+"/preprocessed_log.txt"
     
    # celllogfilepath = outputfolder+"/radmonized_preprocessed_log.txt"
    # Communication_bonesamplefile = outputfolder+"/Communication_radmonized_preprocessed_log.txt"


    
    # outputfolder= './Dataset/Bone/Mouse bone study1/RAW_unlog/'
    
    # # outputfolder= './Dataset/Bone/testsample/'

    # celllogfilepath = outputfolder+"/preprocessed_log.txt"

    # Communication_bonesamplefile = outputfolder+"/Communication_preprocessed_log_v2.txt"



    numthreads=50
    
    l=[]
    

    routenamelist = list(CRoutes.keys())
    samplenumber = 1
    with open(Communication_bonesamplefile,'w') as Communication_bonesample: 
        Header = "Cellname"
        for i in range(len(routenamelist)):
            Header+='\t'+ routenamelist[i] +'-CLR'+'\t'+routenamelist[i] +'-CRR' 
        Communication_bonesample.write(Header+'\n')
        with open(celllogfilepath,'r') as Celldataset:
            # load first row as genename
            istitle = True
            needgeneindex = []
            isstop = False
            # genedict = manager.dict()
            lines = Celldataset.readlines()
            genelist = lines[0].strip().split('\t')
            realgenelist = []
            realgeneorder = []
            for i in range(len(genelist)):
                gene  = genelist[i].upper()
                if gene in RouterelatedGene and not gene in realgenelist:
                    realgenelist.append(gene)
                    realgeneorder.append(i)
            i = 1
            while True:
                writelineslist=manager.list()
                # writelineslist=[]
                cellnamelist =[]
                for j in range(numthreads):
                    if i+j<=len(lines[1:]):
                        linedata = lines[i+j].strip().split('\t')
                        # if linedata[0] in celltypedsamples.keys():
                        #     print(celltypedsamples[linedata[0]])
                        cellnamelist.append(linedata[0])
                        del(linedata[0])
                        samplelist = []
                        for z in realgeneorder:
                            samplelist.append(float(linedata[z]))
                # for CRoutesname,CRoute in CRoutes.items(): 
                        # getRouteline_multiple(writelineslist,samplelist,realgenelist, CRoutes,0.5 ,0.5,i+j,len(lines[1:]))
                        p = Process(target=getRouteline_multiple, args=(writelineslist,samplelist,realgenelist, CRoutes,0.75 ,0.75,i+j,len(lines[1:])))
                        p.start()
                        l.append(p)                  
                    
                for p in l:
                    p.join()

                gc.collect()

                i+=numthreads

                for z in range(len(cellnamelist)):
                    samplenumber+=1
                    print("writing!")
                    print(str(samplenumber))    
                    pline = cellnamelist[z] +"\t"+"\t".join(writelineslist[z])+'\n'
                    Communication_bonesample.write(pline)

                if i>=len(lines[1:]):
                    break
                

        



    os.system("shutdown -s -t  1")
































