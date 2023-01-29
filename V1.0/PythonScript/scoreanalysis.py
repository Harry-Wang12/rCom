




from pyvis.network import Network
import GeneExpressionOperate
import BioknowledgeOperate
import Communication_multiple
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

# def findtopnroute(Probabilityfile,routedetailfile,randonpobabilityfile):

def findthecelltyperesult(celltypeProbailityfile):

    celltyperoutetypecommunication,routelist = GeneExpressionOperate.loadfilewithname(celltypeProbailityfile)

    # start analysis the celltype route score
    # find route wise max
    # find celltype wise max
    # find the celltye wise max -ligand and receptor get violin plot
    # form a network.
    # Get p-value
    
    results ={
        "route wise max":{},
        "celltype wise max":{}
    }

    for route in routelist:
        if route not in results["route wise max"].keys():
            results["route wise max"][route] = {
                "cellpairs":"",
                "value":0
            }

    for celltype,routevalues in celltyperoutetypecommunication.items(): 
        if celltype not in results["celltype wise max"].keys():
            results["celltype wise max"][celltype] = {
                "routename":"",                
                "value":0
            }

            results["celltype wise max"][celltype] ={}


        for routename,routevalue in routevalues.items():
            if len(results["route wise max"][routename]["cellpairs"])==0:
                results["route wise max"][routename]["value"]=routevalue
                results["route wise max"][routename]["cellpairs"] = celltype
            else:
                if routevalue>results["route wise max"][routename]["value"]:
                    results["route wise max"][routename]["value"] = routevalue
                    results["route wise max"][routename]["cellpairs"] = celltype


            results["celltype wise max"][celltype][routename]=routevalue
            
            # else:
            #     if routevalue>results["celltype wise max"][celltype]["value"]:
            #         results["celltype wise max"][celltype]["routename"]=routename
            #         results["celltype wise max"][celltype]["value"] = routevalue

 
    return results


def violinplotforligands(results,rawcellfile,routefile,labelfilepath,sampenamecol,labelnamecol):

    rpacFile = "./DataBaseCombination/KEGG-Pathways/Kegg-Routes.txt"
    pathways_folder= "./DataBaseCombination/KEGG-Pathways/Raw/"
    rpacRoutes,pathwaygraphmap = BioknowledgeOperate.splitrpacRoutes(rpacFile,pathways_folder,"")
    CRoutes,RouterelatedGene = Communication_multiple.loadCommunicationRoutes(routefile,pathwaygraphmap)
    

    # fig,axs = plt.subplots(nrows=subplotnrow,ncols=subplotncol, figsize =(10,6) )


    rawcellsample,genelist = GeneExpressionOperate.loadfilewithname(rawcellfile,transpose=True)
    celllabel = GeneExpressionOperate.loadlabelfile(labelfilepath,sampenamecol,labelnamecol,"\t")


    resultsdf = {}
    for SRpairs,information in results.items():
        routeid = information['routename']
        selectedligands = CRoutes[routeid]['CLR']['Ligand']
        selectedreceptors =  CRoutes[routeid]['CRR']['Receptor']

        scetorname,receivername =  SRpairs.split("_")[0],SRpairs.split("_")[1]
        violinplotdataframe={
            "cellgroup":[],
            "value":[],
            "ligand_receptor":[],
            # "ligandname":[],
            # "recepotorname":[]            
        }

        # for ligand in selectedligands:
        #     if ligand in genelist:
        #         violinplotdataframe[ligand+'-ligand']=[]
        
        # for receptor in selectedreceptors:
        #     if receptor in genelist:
        #         violinplotdataframe[receptor+'-receptor']=[]


        # make first row
        for i in celllabel[scetorname]:
            
            for ligand in selectedligands:
                if ligand in genelist:
                    if i in rawcellsample.keys():
                        violinplotdataframe["cellgroup"].append(scetorname)
                        violinplotdataframe["value"].append(rawcellsample[i][ligand])
                        violinplotdataframe["ligand_receptor"].append("Ligand")

            for receptor in selectedreceptors:
                if receptor in genelist :
                    if i in rawcellsample.keys():
                        violinplotdataframe["cellgroup"].append(scetorname)
                        violinplotdataframe["value"].append(rawcellsample[i][receptor])
                        violinplotdataframe["ligand_receptor"].append("Receptor")    
                
                
        for i in celllabel[receivername]:

            for ligand in selectedligands:
                if ligand in genelist:
                    if i in rawcellsample.keys():
                        violinplotdataframe["cellgroup"].append(receivername)
                        violinplotdataframe["value"].append(rawcellsample[i][ligand])
                        violinplotdataframe["ligand_receptor"].append("Ligand")
            # violinplotdataframe["cellgroup"].append()
            for receptor in selectedreceptors:
                if receptor in genelist :
                    if i in rawcellsample.keys():
                        violinplotdataframe["cellgroup"].append(receivername)
                        violinplotdataframe["value"].append(rawcellsample[i][receptor])
                        violinplotdataframe["ligand_receptor"].append("Receptor")

        resultsdf[SRpairs] = pd.DataFrame.from_dict(violinplotdataframe)

    return resultsdf
        # save file

# def 
 
def generatenetwork(results):
    networkdf = {
        "source":[],
        "target":[],
        "weight":[]
    }

    net = Network(width = "1000px",height ="1000px",notebook=True,directed=True)
    nodenamelist = []
    for SRpairs,information in results.items():

        if not information['value']==1:
        
            routeid = information['routename']
            routevalue = information['value']
            scetorname,receivername =  SRpairs.split("_")[0],SRpairs.split("_")[1]

            if not scetorname in nodenamelist:
                nodenamelist.append(scetorname)
                net.add_node(n_id= nodenamelist.index(scetorname),label=scetorname)
            
            if not receivername in nodenamelist:
                nodenamelist.append(receivername)
                net.add_node(n_id= nodenamelist.index(receivername),label=receivername)

            net.add_edge(source=nodenamelist.index(scetorname),to=nodenamelist.index(receivername),value = routevalue,title = routeid)

        # networkdf["source"].append(scetorname)
        # networkdf["target"].append(receivername)
        # networkdf["weight"].append(routeid+str(routevalue))

    
    # networkdf = pd.DataFrame.from_dict(networkdf)
    # G = nx.from_pandas_edgelist(
    #     networkdf,
    #     source = "source",
    #     target = "target",
    #     edge_attr="weight"
    # )

    # net.from_nx(G)
    net.show_buttons()
    net.show("example.html")

    
def getresult(results,outputfile,routefile,top = 10):

    rpacFile = "./DataBaseCombination/KEGG-Pathways/Kegg-Routes.txt"
    pathways_folder= "./DataBaseCombination/KEGG-Pathways/Raw/"
    rpacRoutes,pathwaygraphmap = BioknowledgeOperate.splitrpacRoutes(rpacFile,pathways_folder,"")
    CRoutes,RouterelatedGene = Communication_multiple.loadCommunicationRoutes(routefile,pathwaygraphmap)




    with open(outputfile,'w') as output:
        output.write("Secetor\tReceptor\tRouteid\tCommunicationValue\tLigand\tReceptor\n") 

        for srpairs,routeinforms in results.items():
            
            sortedrouteinform = {k: v for k, v in sorted(routeinforms.items(), key=lambda item: item[1],reverse=True)}
            
            scetorname,receivername =  srpairs.split("_")[0],srpairs.split("_")[1]
            if not scetorname==receivername:
                selectedlrpairs=[]

                for routename in list(sortedrouteinform.keys()):
                    selectedligands = CRoutes[routename]['CLR']['Ligand']
                    selectedreceptors =  CRoutes[routename]['CRR']['Receptor']
                    if not [selectedligands,selectedreceptors] in selectedlrpairs and len(selectedlrpairs)<=top:
                        selectedlrpairs.append([selectedligands,selectedreceptors])
                        output.write(scetorname+"\t"+receivername+"\t"+routename+"\t"+str(sortedrouteinform[routename])+"\t"+','.join(selectedligands)+"\t"+','.join(selectedreceptors)+"\n") 


if __name__=="__main__":

    # outputfolder= './Dataset/Bone/GSE152285/RAW/'
    # celltypeProbailityfile = outputfolder+"Communication_probability_GSE152285_preprocessed_log_v3.txt"

    # outputfolder= './Dataset/Bone/Mouse bone study1/Raw_processed/'

    # celltypeProbailityfile = outputfolder+"Communication_probability_preprocessed_log_v5.txt"


    outputfolder= './Dataset/Bone/GSE152285/RAW_test6/'
    # celllogfilepath = outputfolder+"/preprocessed.txt"
    # tempfile = outputfolder+"Communication_preprocessed_log_testcorrect_tmp.txt" 
    # Communication_bonesamplefile = outputfolder+'/Communication_preprocessed.txt'

    celltypeProbailityfile = outputfolder+"/Communication_probability_preprocessed_log_v5.txt"


    # labelfilepath = "./Dataset/Bone/GSE152285/meta_annotation_revised.txt"
    # sampenamecol = 0
    # labelnamecol = 1

    labelfilepath = "./Dataset/Bone/Mouse bone study1/mc_compliant_metadata_example.tsv"
    sampenamecol = 0
    labelnamecol = -1

    communicationfilepath = "./DataBaseCombination/RouteStickedResult_V6.txt"

    # rawcellfile = outputfolder+"GSE152285_preprocessed_log.txt"

    CelltypecommunicationResults = findthecelltyperesult(celltypeProbailityfile)


    outputfile = outputfolder+"summary_1.txt"

    getresult(CelltypecommunicationResults["celltype wise max"],outputfile,communicationfilepath,top = 1)

    # with open(outputfile,'w') as output:
    #     output.write("Secetor\tReceptor\tRouteid\tCommunicationValue\n") 

    #     for srpairs,routeinform in CelltypecommunicationResults["celltype wise max"].items():
            
    #         scetorname,receivername =  srpairs.split("_")[0],srpairs.split("_")[1]

    #         output.write(scetorname+"\t"+receivername+"\t"+routeinform['routename']+"\t"+str(routeinform['value'])+"\n") 




    # resultdf = violinplotforligands(CelltypecommunicationResults["celltype wise max"],rawcellfile,communicationfilepath,labelfilepath,sampenamecol,labelnamecol)





    # import seaborn as sns
    # import matplotlib.pyplot as plt

    

    # ax = sns.violinplot(x="cellgroup", y="value",hue="ligand_receptor", data=resultdf["MONOCYTES_SMC"])
    # plt.show()

    # for SRpairs,information in resultdf.items():
    #     print(SRpairs)
    #     plt.figure()
    #     sns.set_theme(style="whitegrid")
    #     tips = sns.load_dataset("tips")

    #     ax = sns.violinplot(x="cellgroup", y="value",hue="ligand_receptor", data=resultdf[SRpairs])

    #     plt.savefig(outputfolder+SRpairs+'.jpg')

   
    # generatenetwork(CelltypecommunicationResults["celltype wise max"])


























