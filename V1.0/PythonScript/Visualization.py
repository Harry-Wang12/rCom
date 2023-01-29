from pyvis.network import Network
import GeneExpressionOperate
import BioknowledgeOperate
import Communication_multiple
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd





def generatenetwork(results,outputfolder):
    networkdf = {
        "source":[],
        "target":[],
        "weight":[]
    }

    net = Network(width = "1000px",height ="1000px",notebook=True,directed=True)
    nodenamelist = []
    for SRpairs,information in results.items():

        if not information['value']==1:
        
            # routeid = information['routename']
            routevalue = information['value']
            scetorname,receivername =  SRpairs.split("_")[0],SRpairs.split("_")[1]

            if not scetorname in nodenamelist:
                nodenamelist.append(scetorname)
                net.add_node(n_id= nodenamelist.index(scetorname),label=scetorname)
            
            if not receivername in nodenamelist:
                nodenamelist.append(receivername)
                net.add_node(n_id= nodenamelist.index(receivername),label=receivername)

            net.add_edge(source=nodenamelist.index(scetorname),to=nodenamelist.index(receivername),value = 1,label = str(routevalue))

        # networkdf["source"].append(scetorname)
        # networkdf["target"].append(receivername)
        # networkdf["weight"].append(routeid+str(routevalue))

    
    # net.from_nx(G)
    net.show_buttons()
    net.show(outputfolder+"example.html")



if __name__=="__main__":
    outputfolder= './Dataset/Bone/GSE152285/'
    # outputfolder= './Dataset/Bone/Mouse bone study1/Raw_processed/'
    # results={
    #     "CC_EC":{
    #         "value":4,
    #     },
    #     "CSM_EC":{
    #         "value":1,
    #     },
    #     "EC_MP":{
    #         "value":2,
    #     },
    #     "EC_OB":{
    #         "value":1,
    #     },
    #     "EC_SMC":{
    #         "value":3,
    #     },
    #     "OB_CC":{
    #         "value":2,
    #     },
    #     "OB_CSM":{
    #         "value":3,
    #     },
    #     "OB_EC":{
    #         "value":6,
    #     },
    #     "OB_MP":{
    #         "value":3,
    #     },
    #     "OB_SBM":{
    #         "value":3,
    #     },
    #     "OB_SMC":{
    #         "value":2,
    #     },
    #     "SBM_EC":{
    #         "value":5,
    #     },
    #     "SMC_EC":{
    #         "value":4,
    #     },

    # }

    results={
            "CC_PL":{
                "value":1,
            },
            "HP_EC":{
                "value":1,
            },
            "HP_PL":{
                "value":3,
            },
            "HP_SMC":{
                "value":1,
            },
            "MC_CC":{
                "value":1,
            },
            "MC_EC":{
                "value":3,
            },
            "MC_OB":{
                "value":3,
            },
            "MC_PH":{
                "value":2,
            },
            "MC_PL":{
                "value":4,
            },
            "MC_SMC":{
                "value":3,
            },
            "OB_HP":{
                "value":1,
            },
            "OB_MC":{
                "value":5,
            },
            "OB_PL":{
                "value":8,
            },
            "PH_MC":{
                "value":1,
            },
            "SMC_MC":{
                "value":9,
            },
            "SMC_PL":{
                "value":11,
            },

    }


    generatenetwork(results,outputfolder)