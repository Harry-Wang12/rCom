
from rPAC import rpac_route as route , rpac_score as score,rpac_pathway_graph as pathway_graph,rpac_util as util
from multiprocessing import Process, Manager
import BioknowledgeOperate
import GeneExpressionOperate
import statistics
import psutil
import gc
import os

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
            Routeid = 'CRid_'+linedata[ResultsHeader.index('id')]
            
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


def getnodeValue(nodetype,nodegenelist,sample,expectation):
    if nodetype == 'gene':
        if expectation>0:
            if nodegenelist[0] in sample.keys():
                if sample[nodegenelist[0]]>0:
                    return sample[nodegenelist[0]]
                else:
                    return 0
            else:
                    return 0

        if expectation<0:
            if nodegenelist[0] in sample.keys():
                if sample[nodegenelist[0]]<0:
                    return abs(sample[nodegenelist[0]])
                else:
                    return 0
            else:
                    return 0
        

    elif nodetype =='And_bundle':
        
        realvaluelist = []
        subweight = 0
        effectivenumber=0
        if expectation>0:
            for genename in nodegenelist:
                if genename in sample.keys(): 
                    effectivenumber+=1                  
                    if sample[genename]>0:

                        realvaluelist.append(sample[genename])
                    if sample[genename]>0:
                        subweight+=1                    
                    elif sample[genename]<0:
                        subweight-=1
            if effectivenumber>0:
                return sum(realvaluelist)/(effectivenumber-subweight+1)
            else:
                return 0

        if expectation<0:
            for genename in nodegenelist:
                if genename in sample.keys(): 
                    effectivenumber+=1                     
                    if sample[genename]<0:
                        realvaluelist.append(sample[genename])
                    if sample[genename]>0:
                        subweight+=1                    
                    elif sample[genename]<0:
                        subweight-=1
            if effectivenumber>0:
                return abs(sum(realvaluelist))/(effectivenumber+subweight+1) 
            else:
                return 0
               
    
    elif nodetype =='OR_bundle':
        genevaluelist = []
        effectivenumber==0       
        for genename in nodegenelist:
            if genename in sample.keys():
                effectivenumber+=1
                genevaluelist.append(sample[genename])

        if expectation>0:     
            if effectivenumber==0:
                return 0
            else:
                if max(genevaluelist)>0:
                    return max(genevaluelist)               
                else:
                    return 0

        if expectation<0:
            for genename in nodegenelist:
                if genename in sample.keys(): 
                    effectivenumber+=1                     
                    if sample[genename]<0:
                        realvaluelist.append(sample[genename])
                    if sample[genename]>0:
                        subweight+=1                    
                    elif sample[genename]<0:
                        subweight-=1
            if effectivenumber>0:
                return abs(sum(realvaluelist))/(len(effectivenumber)+subweight+1) 
            else:
                return 0 


def CLRcalculation(sample,CLR,weightratio = 0.75):
    # CLR = {        
    #     'Up':[],
    #     'Ligand':[],
    #     'Expectation':[],
    #     'Types':[]
    # }
    # the value suppose get value from 0-1
    Ligandvalue = getnodeValue(CLR['Types'][-1],CLR['Ligand'],sample,1)
    Upvalue = 0
    for i in range(len(CLR['Up'])):
        Upvalue+=  getnodeValue(CLR['Types'][i],CLR['Up'][i],sample,CLR['Expectation'][i])

    return weightratio*Ligandvalue+Upvalue*(1-weightratio)

def CRRcalculation(sample,CRR,weightratio = 0.75):
    # CRR = {
    #     'Receptor':[],        
    #     'Down':[],        
    #     'Expectation':[],
    #     'Types':[]        
    # }

    # the value suppose get value from 0-1

    Receptorvalue = getnodeValue(CRR['Types'][0],CRR['Receptor'],sample,1)
    Downvalue = 0
    for i in range(len(CRR['Down'])):
        Downvalue+= getnodeValue(CRR['Types'][i+1],CRR['Down'][i],sample,CRR['Expectation'][i+1])
        
    return weightratio*Receptorvalue+Downvalue*(1-weightratio)


def calculateCLRandCRR(CommunicationRoute,sample,CLRweight = 0.75,CRRweight = 0.75 ):
    
    CLRscore = CLRcalculation(sample,CommunicationRoute['CLR'],CLRweight)
    CRRscore = CRRcalculation(sample,CommunicationRoute['CRR'],CRRweight)

    return CLRscore,CRRscore

def getcellcommunicationVector(cellsample, CommunicationRoutes,CLRweight ,CRRweight):
    newsample = {}
    for routeid,CommunicationRoute in CommunicationRoutes.items():
        
        # if(routeid == "CommunicationRoueid_157"):
        #     print("11")
        CLRscore,CRRscore = calculateCLRandCRR(CommunicationRoute,cellsample,CLRweight ,CRRweight)
        newsample[routeid+'-CLR'] = CLRscore
        newsample[routeid+'-CRR'] = CRRscore
    # newsamples[samplename] = newsample
    # print("Done with "+ samplename)
    return newsample


def loaddataset_multiple( inputfilepath, relatedgene,sep = '\t',transpose=False):

    genelist = []
    allnamesamples={}
    with open(inputfilepath,'r') as allsamplefile:
        istitle=True
        if transpose:
            genelist = [] 
            for line in allsamplefile.readlines():
                linedata = line.strip().split(sep)
                if istitle:
                    del(linedata[0])
                    for genename in linedata:
                        if genename in relatedgene:
                            genelist.append(genename.upper())                        
                    istitle=False
                else:
                    sampledetails ={} 
                    del(linedata[0])
                    for i in range(len(genelist)):                        
                        sampledetails[genelist[i]] = float(linedata[i])
                    allnamesamples.append(sampledetails)                    
        else:
            samplelist=[]
            for line in allsamplefile.readlines():
                linedata = line.strip().split(sep)
                if istitle:
                    del(linedata[0])
                    for sample in linedata:
                        allnamesamples[sample]=[]
                        samplelist.append(sample)
                    istitle=False
                else:                    
                    genename = linedata[0].upper()
                    if genename not in genelist:
                        if genename in relatedgene:
                            genelist.append(genename)
                    del(linedata[0])
                    for i in range(len(samplelist)):
                        allnamesamples[samplelist[i]].append(float(linedata[i])) 
                        
    # print("dup")
    # mgenelist = Mmanager.list()
    # mBonecells= Mmanager.dict()

    
    # for cellname,cell in allnamesamples.items():
    #     mBonecells[cellname] = cell

    # for genename in genelist:
    #     mgenelist.append(genename)

    # print("done")


    # return mBonecells,mgenelist
    return allnamesamples,genelist





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
                    return abs(sample[Locindex])
                else:
                    return 0
            else:
                    return 0
        

    elif nodetype =='And_bundle':
        
        realvaluelist = []
        subweight = 0
        effectivenumber=0
        if expectation>0:
            for genename in nodegenelist:
                if genename in genelist: 
                    effectivenumber+=1  
                    Locindex = genelist.index(genename)                
                    if sample[Locindex]>0:
                        realvaluelist.append(sample[Locindex])
                    if sample[Locindex]>0:
                        subweight+=1                    
                    elif sample[Locindex]<0:
                        subweight-=1
            if effectivenumber>0:
                return sum(realvaluelist)/(effectivenumber-subweight+1)
            else:
                return 0

        if expectation<0:
            for genename in nodegenelist:
                if genename in genelist: 
                    effectivenumber+=1  
                    Locindex = genelist.index(genename)                   
                    if sample[Locindex]<0:
                        realvaluelist.append(sample[Locindex])
                    if sample[Locindex]>0:
                        subweight+=1                    
                    elif sample[Locindex]<0:
                        subweight-=1
            if effectivenumber>0:
                return abs(sum(realvaluelist))/(effectivenumber+subweight+1) 
            else:
                return 0
               
    
    elif nodetype =='OR_bundle':
        realvaluelist = []
        subweight = 0
        genevaluelist = []
        effectivenumber=0       
        for genename in nodegenelist:
            if genename in genelist:
                Locindex = genelist.index(genename)      
                effectivenumber+=1
                genevaluelist.append(sample[Locindex])

        if expectation>0:     
            if effectivenumber==0:
                return 0
            else:
                if max(genevaluelist)>0:
                    return max(genevaluelist)               
                else:
                    return 0

        if expectation<0:
            for genename in nodegenelist:
                if genename in genelist: 
                    Locindex = genelist.index(genename) 
                    effectivenumber+=1                     
                    if sample[Locindex]<0:
                        realvaluelist.append(sample[Locindex])
                    if sample[Locindex]>0:
                        subweight+=1                    
                    elif sample[Locindex]<0:
                        subweight-=1
            if effectivenumber>0:
                return abs(sum(realvaluelist))/(effectivenumber+subweight+1) 
            else:
                return 0 
    







def getRouteSampleline_multiple(resultlist, cellsamples,genelist, Routename ,CommunicationRoute,CLRweight ,CRRweight):
    print(Routename)
    
        
    # if(routeid == "CommunicationRoueid_157"):
    #     print("11")

    CLRline = Routename+'-CLR'
    CRRline = Routename+'-CRR'


    for samplename,samplevaluelist in cellsamples.items():
        # CLR
        # 
        # 
    #    getnodeValue_mulitple(nodetype,nodegenelist,sample,expectation,genelist)
        Ligandvalue = getnodeValue_mulitple(CommunicationRoute['CLR']['Types'][-1],CommunicationRoute['CLR']['Ligand'],samplevaluelist,1,genelist)
        Upvalue = 0
        for i in range(len(CommunicationRoute['CLR']['Up'])):
            Upvalue+=  getnodeValue_mulitple(CommunicationRoute['CLR']['Types'][i],CommunicationRoute['CLR']['Up'][i],samplevaluelist,CommunicationRoute['CLR']['Expectation'][i],genelist)

        CLRline+="\t"+str( CLRweight*Ligandvalue+Upvalue*(1-CLRweight))

        # CLRscore =CLRcalculation(cellsample,CommunicationRoute['CLR'],CLRweight)
        
        # CRR

        Receptorvalue = getnodeValue_mulitple(CommunicationRoute['CRR']['Types'][0],CommunicationRoute['CRR']['Receptor'],samplevaluelist,1,genelist)
        Downvalue = 0
        for i in range(len(CommunicationRoute['CRR']['Down'])):
            Downvalue+= getnodeValue_mulitple(CommunicationRoute['CRR']['Types'][i+1],CommunicationRoute['CRR']['Down'][i],samplevaluelist,CommunicationRoute['CRR']['Expectation'][i+1],genelist)
            
        CRRline+="\t"+str(CRRweight*Receptorvalue+Downvalue*(1-CRRweight))


        # CRRscore = CRRcalculation(cellsample,CommunicationRoute['CRR'],CRRweight)

    resultlist.append(CLRline+'\n')
    resultlist.append(CRRline+'\n')
    print("Done "+str(len(resultlist)))
    # newsamples[samplename] = newsample
    # print("Done with "+ samplename)
    # return newsample



def getRouteline_multiple(resultlist, sample,genelist, CommunicationRoutes,CLRweight ,CRRweight,Cindex,totalnumber):
    
    print(str(Cindex)+'/'+str(totalnumber))    
    samplelist = []


    for routename,CommunicationRoute in CommunicationRoutes.items():
        # CLR
        # 
        # 
    #    getnodeValue_mulitple(nodetype,nodegenelist,sample,expectation,genelist)
        Ligandvalue = getnodeValue_mulitple(CommunicationRoute['CLR']['Types'][-1],CommunicationRoute['CLR']['Ligand'],sample,1,genelist)
        Upvalue = 0
        for i in range(len(CommunicationRoute['CLR']['Up'])):
            Upvalue+=  getnodeValue_mulitple(CommunicationRoute['CLR']['Types'][i],CommunicationRoute['CLR']['Up'][i],sample,CommunicationRoute['CLR']['Expectation'][i],genelist)

        samplelist.append (str(round(CLRweight*Ligandvalue+Upvalue*(1-CLRweight),2)))

        # CLRscore =CLRcalculation(cellsample,CommunicationRoute['CLR'],CLRweight)
        
        # CRR

        Receptorvalue = getnodeValue_mulitple(CommunicationRoute['CRR']['Types'][0],CommunicationRoute['CRR']['Receptor'],sample,1,genelist)
        Downvalue = 0
        for i in range(len(CommunicationRoute['CRR']['Down'])):
            Downvalue+= getnodeValue_mulitple(CommunicationRoute['CRR']['Types'][i+1],CommunicationRoute['CRR']['Down'][i],sample,CommunicationRoute['CRR']['Expectation'][i+1],genelist)
            
        samplelist.append(str(round(CRRweight*Receptorvalue+Downvalue*(1-CRRweight),2)))


        # CRRscore = CRRcalculation(cellsample,CommunicationRoute['CRR'],CRRweight)

    resultlist.append(samplelist)
    # resultlist.append(CRRline+'\n')
    print("Done "+str(Cindex)+'/'+str(totalnumber))
    # newsamples[samplename] = newsample
    # print("Done with "+ samplename)
    # return newsample



def CalculateCoummunicationMatrixtolist(CommunicationRoutes,samples,CLRweight = 0.75,CRRweight = 0.75 ):

    outputlist = []
    header = "Comunication_ID\t"+"\t".join(list(samples.keys()))+'\n'
    outputlist.append(header)

    for routeid,CommunicationRoute in CommunicationRoutes.items():
        CLRline = routeid+'-CLR'
        CRRline = routeid+'-CRR'
        for samplename,cellsample in samples.items():
            CLRscore,CRRscore = calculateCLRandCRR(CommunicationRoute,cellsample,CLRweight ,CRRweight)
            CLRline+='\t'+str(CLRscore)
            CRRline+='\t'+str(CRRscore)
        
        outputlist.append(CLRline+'\n')
        outputlist.append(CRRline+'\n')
        print(routeid)
        # if routeid =='CommunicationRoueid_10':
        #     return outputlist
    return outputlist



# def CalculateCoummunicationMatrixtofile(CommunicationRoutes,samples,CLRweight = 0.75,CRRweight = 0.75, outputfile="" ) :
#     # use it when memory is low
#     newsamples={}
#     with Manager() as manager:
#         newsamples=manager.dict()
#         l = []
#         for samplename,sample in samples.items():
#             # change to multiple process
#             p = Process(target=getcellcommunicationVector, args=(newsamples,samplename,sample, CommunicationRoutes,CLRweight ,CRRweight,))
#             p.start()
#             l.append(p)
#             # p.join()
#         for p in l:
#             p.join()
        
    # for samplename,sample in samples.items():
    #     newsamples[samplename] = getcellcommunicationVector(sample, CommunicationRoutes,CLRweight ,CRRweight)


    # with open(outputfile,'w') as output:
    #     header = "Comunication_ID\t"+"\t".join(list(samples.keys()))+'\n'
    #     output.write(header)
    #
    #     for routeid,CommunicationRoute in CommunicationRoutes.items():
    #         CLRline = routeid+'-CLR'
    #         CRRline = routeid+'-CRR'
    #         for samplename,cellsample in samples.items():
    #             CLRscore,CRRscore = calculateCLRandCRR(CommunicationRoute,cellsample,CLRweight ,CRRweight)
    #             CLRline+='\t'+str(CLRscore)
    #             CRRline+='\t'+str(CRRscore)
    #
    #         output.write(CLRline+'\n')
    #         output.write(CRRline+'\n')
    #         print(routeid)


    # return newsamples







if __name__=='__main__':



    manager = Manager() 
    ChipSeqfilepath = "./DataBaseCombination/ChipSeqLibrary/gene_attribute_matrix.txt"

    ChipSeqDB= BioknowledgeOperate.loadChipSeqfile(ChipSeqfilepath)
    print('RAM memory % used:', psutil.virtual_memory()[2])

    CellTalkfilepath = "./DataBaseCombination/CellTalkLRDB/human_lr_pair.txt"

    CellTalkDB = BioknowledgeOperate.loadCellTalkDB(CellTalkfilepath)
    print('RAM memory % used:', psutil.virtual_memory()[2])

    CellChatfilepath = "./DataBaseCombination/CellChatLRDB/Human_complex_input_CellChatDB.csv"

    CellChatDB = BioknowledgeOperate.loadCellChatDB(CellChatfilepath)
    print('RAM memory % used:', psutil.virtual_memory()[2])

    outputfile = "./DataBaseCombination/RouteStickedResult_V5.txt"

    rpacFile = "./DataBaseCombination/KEGG-Pathways/Kegg-Routes.txt"

    pathways_folder= "./DataBaseCombination/KEGG-Pathways/Raw/"

    rpacRoutes,pathwaygraphmap = BioknowledgeOperate.splitrpacRoutes(rpacFile,pathways_folder,"")

    CRoutes,RouterelatedGene = loadCommunicationRoutes(outputfile,pathwaygraphmap)


    # foloderlist = ['SMC','Platelets','MSC','Proliferating_Hem','Hematopoietic','Osteoblast','Monocytes']
# 
    # foloderlist = ['Monocytes',]

    

    # for folder in foloderlist:
    # print(folder)
    outputfolder= './Dataset/Bone/Normal/Normal3/produced'
    bonesamplefile = outputfolder+"/cart_normal1_preprocessed_log.txt"

    Communication_bonesamplefile = outputfolder+"/Communication_cart_normal1_preprocessed_log_0.5.txt"

    iswrite = True

    # socketPairs= BioknowledgeOperate.analysisSocketPairs(ChipSeqDB, CellTalkDB,CellChatDB)

    # BioknowledgeOperate.generateCellCommunicationDB_V1(rpacRoutes,pathwaygraphmap,socketPairs,outputfile,iswrite)

    
    print('RAM memory % used:', psutil.virtual_memory()[2])

    Bonecells,genelist = loaddataset_multiple(bonesamplefile,RouterelatedGene,'\t')
    print('RAM memory % used:', psutil.virtual_memory()[2])
    
    
    # mBonecells = manager.dict()
    # for cellname,cell in Bonecells.items():
    #     mBonecells[cellname] = cell

    # del(Bonecells)
    # gc.collect()
    print('RAM memory % used:', psutil.virtual_memory()[2])
    # del(pathwaygraphmap)
    
    writelineslist=manager.list()
    l = []
    
    # writelineslist.append(header)
    numthreads= 30
    Cellnamelist = list(Bonecells.keys())
    i=0
    # getRouteline_multiple(writelineslist, Bonecells[Cellnamelist[1]],genelist, CRoutes,0.75 ,0.75,1,10)

###########################################################
    while True:

        for j in range(numthreads):
            if i+j<len(Cellnamelist):
        # for CRoutesname,CRoute in CRoutes.items(): 
                p = Process(target=getRouteline_multiple, args=(writelineslist, Bonecells[Cellnamelist[i+j]],genelist, CRoutes,0.5 ,0.5,i+j,len(Cellnamelist)))
                p.start()
                l.append(p)                  
                del(Bonecells[Cellnamelist[i+j]])
               

        for p in l:
            p.join()
        gc.collect()
        i+=numthreads
        if i>=len(Cellnamelist):
            break
###########################################################
# 
#       
    # i = 1
    # for cellname in Cellnamelist:
    #     cellsample= Bonecells[cellname]
    #     p = Process(target=getRouteline_multiple, args=(writelineslist, cellsample,genelist, CRoutes,0.5 ,0.5,i,len(Cellnamelist)))
    #     i+=1
    #     p.start()    
    #     l.append(p)
    
    # for p in l:
    #     p.join()

    print('RAM memory % used:', psutil.virtual_memory()[2])
    del(Bonecells)
    gc.collect()
    print('RAM memory % used:', psutil.virtual_memory()[2])





    # # writeline=[]
    # # print("writing")
    # # X = 0
    # # for sample in writelineslist:
    # #     X+=1
    # #     print(X)
    # #     writeline.append(list(sample))
    
    
    # print('RAM memory % used:', psutil.virtual_memory()[2])
    # del(writelineslist)
    # gc.collect()
    # print('RAM memory % used:', psutil.virtual_memory()[2])
    
    # header = "Comunication_ID\t"+"\t".join(list(Bonecells.keys()))+'\n'

    routenamelist = list(CRoutes.keys())

    with open(Communication_bonesamplefile,'w') as Communication_bonesample: 
        Header = "Cellname"
        for i in range(len(routenamelist)):
            Header+='\t'+ routenamelist[i] +'-CLR'+'\t'+routenamelist[i] +'-CRR' 
        Communication_bonesample.write(Header+'\n')

        
        for i in range(len(Cellnamelist)):
            print(str(i)+'/'+str(len(Cellnamelist)))    
            pline = Cellnamelist[i] +"\t"+"\t".join(writelineslist[i])+'\n'
            Communication_bonesample.write(pline)





        # with open(Communication_bonesamplefile,'w') as Communication_bonesample: 
        #     Communication_bonesample.write(header) 

        #     for i in range(len(routenamelist)):
        #         print( routenamelist[i])
        #         CLRline = routenamelist[i] +'-CLR' 
        #         CRRline = routenamelist[i] +'-CRR' 
        #         CLRindex = i*2
        #         CRRindex = i*2+1
                
        #         for sample in writeline:
        #             CLRline+='\t'+str(sample[CLRindex])
        #         Communication_bonesample.write(CLRline+'\n') 

        #         for sample in writeline:
        #             CRRline+='\t'+str(sample[CRRindex])
        #         Communication_bonesample.write(CRRline+'\n')  


    # os.system("shutdown -s -t  1")


    # writeallline = CalculateCoummunicationMatrixtolist(CRoutes,Bonecells,CLRweight = 0.75,CRRweight = 0.75 )
    
    # with open(Communication_bonesamplefile,'w') as Communication_bonesample: 
    #     for line in writeallline:
    #         Communication_bonesample.write(line)  

    # os.system("shutdown -s -t  1")
    # index=0
    # with Manager() as manager:
    #     newsamples=manager.dict()
    #     l = []
    #     for samplename,sample in Bonecells.items(): 
    #         index+=1
    #         # change to multiple process    
    #         p = Process(target=getcellcommunicationVector, args=(newsamples,samplename,sample, CRoutes,0.75 ,0.75,index))
    #         p.start()
    #     #     l.append(p)
    #         # p.join() 
    #     for p in l:
    #         p.join()

    # GeneExpressionOperate.writenamedsampletofile(newsamples,Communication_bonesamplefile)
























