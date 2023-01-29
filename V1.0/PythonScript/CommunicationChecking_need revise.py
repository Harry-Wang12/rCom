
from socket import socketpair
from zmq import Socket
import GeneExpressionOperate

def AddingRouteinformationToSamples(GeneExpressionDatasetfile,RouteScoringFile):
    
    '''
        input 
            GeneExpressionDatasetfile,RouteScoringFile
        output
            Samples:
                genevalue
                RouteScoringFile
    '''
    Genesamples,genelist = GeneExpressionOperate.loadfilewithname(GeneExpressionDatasetfile)
    Routesamples,routelist = GeneExpressionOperate.loadfilewithname(RouteScoringFile)


    CombinedCells={}

    for Cellid, Geneinformation in Genesamples.items():
        CombinedCells[Cellid]= Geneinformation
        for routeid,RouteScore in Routesamples[Cellid].items():
            CombinedCells[Cellid][routeid] = RouteScore



    return CombinedCells


def IdentifiedLRCells(LRBDfile,CombinedCells,outputfolder,routeThreshold=0.1, geneThreshold=0.1):
    ResultsHeader = [
        'id',
        'C2-R2-ID',
        'C2-R2-PW',
        'C2-R2-TYPE',
        'C2-R2-ids',
        'C2-R2-genes',
        'C2-R2-RRecetors-ids',
        'C2-R2-RRecetors-genes',
        'C2-R2-RLigands-ids',
        'C2-R2-RLigands-genes',
        'SocketPairs_Receptor',
        'SocketPairs_Ligand',
        'SocektPairs_TF',
        'C1-R1-TFs-ids',
        'C1-R1-TFs-genes',
        'C1-R1-ids',
        'C1-R1-genes',
        'C1-R1-PW',
        'C1-R1-ID',
        'C1-R1-TYPE',
        'combineType']
    
    istitle = False
    CommunicationResult = {}
    with open(LRBDfile,'r') as LRDB:
        for line in LRDB.readlines():
            if istitle:
                istitle =False
            else:
                linedata = line.strip().split('\t')
                outputfile = outputfolder+'Combination_'+linedata[ResultsHeader.index('id')]+'.txt'
                print(outputfile)
                CommunicationResult['Combination_'+linedata[ResultsHeader.index('id')]] = {
                    'Receivers':[],
                    'Secretors':[]
                }
                CombinationType = linedata[ResultsHeader.index('combineType')]
                with open(outputfile,'w') as outf:
                                
                    if CombinationType =='4':

                        C2_R2_ID = linedata[ResultsHeader.index('C2-R2-ID')].upper()

                        socketpairTF = linedata[ResultsHeader.index('SocektPairs_TF')]
                        socketpairLigand = linedata[ResultsHeader.index('SocketPairs_Ligand')]
                       
                        for sampleid,sampleValues in CombinedCells.items():
                            if sampleValues[C2_R2_ID] >=routeThreshold:
                                CommunicationResult['Combination_'+linedata[ResultsHeader.index('id')]]['Receivers'].append(sampleid)
                                outf.write("Receivers:"+sampleid+'\n')
                            if socketpairLigand in CombinedCells.keys():
                                if socketpairTF in CombinedCells.keys():
                                    if sampleValues[socketpairTF] >=geneThreshold and sampleValues[socketpairLigand] >=geneThreshold:
                                        CommunicationResult['Combination_'+linedata[ResultsHeader.index('id')]]['Secretors'].append(sampleid)
                                        outf.write("Secretors:"+sampleid+'\n')
                                elif sampleValues[socketpairLigand] >=geneThreshold:
                                    CommunicationResult['Combination_'+linedata[ResultsHeader.index('id')]]['Secretors'].append(sampleid)
                                    outf.write("Secretors:"+sampleid+'\n')
                    if CombinationType =='5':

                        C2_R2_ID = linedata[ResultsHeader.index('C2-R2-ID')].upper()

                        socketpairTF = linedata[ResultsHeader.index('SocektPairs_TF')]
                        socketpairLigand = linedata[ResultsHeader.index('SocketPairs_Ligand')]

                        for sampleid,sampleValues in CombinedCells.items():
                            if sampleValues[C2_R2_ID] >=routeThreshold:
                                CommunicationResult['Combination_'+linedata[ResultsHeader.index('id')]]['Receivers'].append(sampleid)
                                outf.write("Receivers:"+sampleid+'\n')
                            if socketpairLigand in CombinedCells.keys():
                                if sampleValues[socketpairLigand] >=geneThreshold:
                                    CommunicationResult['Combination_'+linedata[ResultsHeader.index('id')]]['Secretors'].append(sampleid)
                                    outf.write("Secretors:"+sampleid+'\n')
                    
                    if CombinationType =='2':
                        C2_R2_ID = linedata[ResultsHeader.index('C2-R2-ID')].upper()
                        C1_R1_ID = linedata[ResultsHeader.index('C1-R1-ID')].upper()
                        socketpairTF = linedata[ResultsHeader.index('SocektPairs_TF')]
                        socketpairLigand = linedata[ResultsHeader.index('SocketPairs_Ligand')]

                        for sampleid,sampleValues in CombinedCells.items():
                            if sampleValues[C2_R2_ID] >=routeThreshold:
                                CommunicationResult['Combination_'+linedata[ResultsHeader.index('id')]]['Receivers'].append(sampleid)
                                outf.write("Receivers:"+sampleid+'\n')
                            if sampleValues[C1_R1_ID] >=routeThreshold:
                                CommunicationResult['Combination_'+linedata[ResultsHeader.index('id')]]['Secretors'].append(sampleid)
                                outf.write("Secretors:"+sampleid+'\n')

                


    return CommunicationResult        

                





































