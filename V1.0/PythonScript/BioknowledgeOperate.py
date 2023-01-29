import os
import xml.etree.ElementTree as ET


from rPAC import rpac_route as route , rpac_score as score,rpac_pathway_graph as pathway_graph,rpac_util as util

# import rPAC.rpac_route as route 
# import rPAC.rpac_score as score
# import rPAC.rpac_pathway_graph as pathway_graph
# import rPAC.rpac_util as util

def loadTrrust(filepath, includeregulation=['Repression','Unknown','Activation']):
    '''
    load Trrust database from Trrust file
    the file have three column 
    TF \t target \t regrulation \n

    includeregulation = ['Repression','Unknown','Activation']

    return edgelist[
        'Causal': tfgenename,
        'Target': targetgenename,
        'Regulation': regulation
    ]
    '''

    edgelist=[]

    with open(filepath,'r') as trrustfile:
        for line in trrustfile.readlines():
            linedata= line.strip().split("\t")
            
            tfgenename = linedata[0].upper()
            targetgenename = linedata[1].upper()
            regulation = linedata[2]
            if regulation in includeregulation:
                edgelist.append({
                    'Causal': tfgenename,
                    'Target': targetgenename,
                    'Regulation': regulation
                })

    return edgelist



    
def getAllCausalGenesforOneTarget(edgelist,Targetgene):
    
    '''
    Return all the causal genes for the target
    input:
        edgelist
        target gene name
    '''
    causalgenes = []
    for edge in edgelist:
        if edge['Target']==Targetgene:
            causalgenes.append(edge['Causal'])
    return causalgenes



def getAllCausalGenesforAllTargets(edgelist):

    '''
    Return all the causal genes for all targets
    input:
        edgelist        
    '''

    edgedict={}
    for edge in edgelist:
        if not edge['Target'] in edgedict.keys():
            edgedict[edge['Target']] = []
        
        edgedict[edge['Target']].append(edge['Causal'])

    return edgedict


def getEdgeRegulation(edgelist, Causalgene, Targetgene):
    
    '''
    Get the specific edge regulation from edgelist
    '''
    for edge in edgelist:
       if edge['Target']==Targetgene and edge['Causal']==Causalgene:
           return  edge['Regulation'] 

def parsarKeggfromFolder(folderpath):
    kegglist = {}
    g = os.walk(folderpath)
    for path, dir_list,file_list in g:
        for file_name in file_list:
            pathwayname,nodes,relation,genepairs = parsarKeggintonodepairs(os.path.join(path,file_name))
            kegglist[pathwayname] = {
                "nodes":nodes,
                "relation":relation,
                "genepairs":genepairs
            }
   
    return kegglist
        
def parsarKeggintonodepairs(kgmlfilepath):

    '''
    Get the genepairs from KGML file
    '''
    # print(kgmlfilepath)
    tree = ET.parse(kgmlfilepath)
    root = tree.getroot()
    pathwayname = root.attrib["title"]
    nodes = {}
    relation = []
    genepairs = []
    for child in root:
        if child.tag == 'entry' and not child.attrib['name'] == 'undefined':
            id = child.attrib['id']
            nodetype = child.attrib['type']
            for subchild in child:
                namelist =subchild.attrib['name'].replace("...","").strip().split(',')
            nodes[id]={
            "nodetype":nodetype,
            "namelist":namelist 
            }
        if child.tag == 'entry' and child.attrib['name'] == 'undefined':
            id = child.attrib['id']
            nodetype = child.attrib['type']
            namelist = []
            for subchild in child:
                if subchild.tag == "component":
                    namelist+=nodes[subchild.attrib['id']]["namelist"]
            nodes[id]={
            "nodetype":nodetype,
            "namelist":namelist 
            }
        if child.tag == 'relation' :
            start = child.attrib['entry1']
            end = child.attrib['entry2']
            if len(list(child)) :
                type = child[0].attrib['name']
            else:
                type = None
            relation.append({
               "start":start,
               "end":end ,
               "type":type
            })
            if start in nodes.keys() and end in nodes.keys() :
                for gene1 in nodes[start]['namelist']:
                    for gene2 in nodes[end]['namelist']:
                        genepairs.append([gene1,gene2,type])
    return pathwayname,nodes,relation,genepairs

def analysisroute(genepairs, routegenes):
    
    routegenepairs = []
    for genepair in genepairs:
        gene1 = genepair[0]
        gene2 = genepair[1]
        if gene1 in routegenes and gene2 in routegenes:
            routegenepairs.append(genepair)
    return routegenepairs






def Formpathwaycomponets(KEGGxmlfolder,KEGGjsonfolder):

    # get route
    ligand_list = ["LIGAND", "BUNDLE_LIGAND", "BUNDLE_OR_LIGAND",
               "BUNDLE_AND_LIGAND", "BUNDLELIGANDAND", "BUNDLELIGANDOR"]
    receptor_list = ["RECEPTOR", "BUNDLE_RECEPTOR",
                    "BUNDLE_OR_RECEPTOR", "BUNDLE_AND_RECEPTOR"]
    tf_list = ["TRANSCRIPTION FACTOR", "BUNDLE_AND_TF", "BUNDLE_OR_TF",
            "TRANSCRIPTIONFACTOR", "DIAMOND", "BUNDLETFCAND", "BUNDLETFCOR"]
    bp_list = ["BIOLOGICAL PROCESS"]
    pathway_type_list = ['PATHWAY']

    route_end_candidate_types = {'P1': {'SOURCE': ligand_list + receptor_list, 'TARGET': tf_list + pathway_type_list},
                             'P2': {'SOURCE': tf_list, 'TARGET': [], 'LEAF_ALLOWED': True}}


    all_routes = route.find_routes_on_pathways(KEGGjsonfolder, route_end_candidate_types, False, "")   
    keggdict = parsarKeggfromFolder(KEGGxmlfolder)

    for pathwayroute in all_routes[1:]:        
        pathwayname = pathwayroute[3]
        routegenes = pathwayroute[5]   
        routeidex = pathwayroute[0]     
        if pathwayname not in keggdict.keys():
            print(pathwayname)
        else:
            if 'routes' not in keggdict[pathwayname].keys():
                keggdict[pathwayname]['routes']={}
            keggdict[pathwayname]['routes'][routeidex]= analysisroute(keggdict[pathwayname]['genepairs'], routegenes)
    
    # get the genes pairs that is not included in the routes
    
    for pathwayname,pathwayinfo in keggdict.items():
        # get genes included in the routes
        selectedroutegenepairs = []
        functionalgroupgenepairs=[]
        if 'routes' not in pathwayinfo.keys():
            functionalgroupgenepairs = pathwayinfo['genepairs']
        else:
            for routename, routegenepairs in pathwayinfo['routes'].items():
                for genepair in routegenepairs:
                    if genepair not in selectedroutegenepairs:
                        selectedroutegenepairs.append(genepair)
            
            for genepair in pathwayinfo['genepairs']:

                if genepair not in selectedroutegenepairs:
                        functionalgroupgenepairs.append(genepair)


        keggdict[pathwayname]['prefunctionalgroup']=functionalgroupgenepairs

    return keggdict

    # print(all_routes[0])
    # find the matched pathway information
    # Calcualted detailed route in this pathway.


def analysisSocketPairs(ChipSeqDB, CellTalkDB,CellChatDB):
    # header_row  = ["ID","TF","Ligand(Target)","Receptor","Source"]
    # result_rows = [header_row]
    # collection = []
    collectiondict = {}
    # index= 1
    for receptor,ligands in CellTalkDB.items():
        if receptor not in collectiondict.keys():
            collectiondict[receptor]={}
        for ligand in ligands:
            if ligand not in collectiondict[receptor].keys():
                collectiondict[receptor][ligand]=[]
            if ligand in ChipSeqDB.keys():
                for tf in ChipSeqDB[ligand]:
                    if tf not in collectiondict[receptor][ligand]:
                        collectiondict[receptor][ligand].append(tf)
    for receptor,ligands in CellChatDB.items():
        if receptor not in collectiondict.keys():
            collectiondict[receptor]={}
        for ligand in ligands:
            if ligand not in collectiondict[receptor].keys():
                collectiondict[receptor][ligand]=[]
            if ligand in ChipSeqDB.keys():
                for tf in ChipSeqDB[ligand]:
                    if tf not in collectiondict[receptor][ligand]:
                        collectiondict[receptor][ligand].append(tf)

    return collectiondict

















    # for tf, targetlist in ChipSeqDB.items():    
    #     # print(tf)    
    #     for targetgene in targetlist:
    #         if targetgene in CellTalkDB.keys():
    #             ligand = targetgene
    #             for receptor in CellTalkDB[ligand]:
    #                 if [[tf,ligand,receptor]] not in collection:
    #                     collection.append([tf,ligand,receptor])
    #                     result_rows.append([str(index),tf,ligand,receptor,"CellTalk"])

    #                     if receptor not in collectiondict.keys():
    #                         collectiondict[receptor]={}
    #                     if ligand not in collectiondict[receptor].keys():
    #                         collectiondict[receptor][ligand]=[]
    #                     if tf not in collectiondict[receptor][ligand]:
    #                         collectiondict[receptor][ligand].append(tf)

    #                     index +=1
    #         if targetgene in CellChatDB.keys():
    #             ligand = targetgene
    #             for receptor in CellChatDB[ligand]:
    #                 if [[tf,ligand,receptor]] not in collection:
    #                     collection.append([tf,ligand,receptor])
    #                     result_rows.append([str(index),tf,ligand,receptor,"CellChat"])
    #                     if receptor not in collectiondict.keys():
    #                         collectiondict[receptor]={}
    #                     if ligand not in collectiondict[receptor].keys():
    #                         collectiondict[receptor][ligand]=[]
    #                     if tf not in collectiondict[receptor][ligand]:
    #                         collectiondict[receptor][ligand].append(tf)
    #                     index +=1
    # if iswrite:
    #     with open(outputfile,'w') as f:
    #         for row in result_rows:
    #             f.write('\t'.join(row)+'\n')
    



# def formfunctionalgroup(genepairs):
    





    # return all_routes




def loadChipSeqfile(filepath):
    TFlist = []
    targetdict = {}

    with open(filepath, 'r') as f:
        linenumber = 0
        # 1st line 
            # [0] = #
            # [1] = #
            # [2] = GeneSym
            # [3:] = genesyms
        # 2nd line not use
        # 3rd line not use
        # 4th: line start 0 means not target 1 means target
        for line in f.readlines():
            linenumber+=1
            
            if linenumber == 1:
                linedata = line.strip().split('\t')
                TFlist = linedata[3:]
            if linenumber >3:
                linedata = line.strip().split('\t')
                targetname = linedata[0]
                belonglist = linedata[3:]
                for i in range(len(belonglist)):
                    if belonglist[i] == "1.000000":
                        if targetname not in targetdict.keys():
                            targetdict[targetname] = []
                        targetdict[targetname].append(TFlist[i])

    return targetdict


def loadCellChatDB(lrfile, ligandcol=1, receptorcol=2):
    lrdict={}
    with open(lrfile,'r') as f:
        title = True
        for line in f.readlines():
            
            if title:
                title = False
            else:
                linedata = line.replace('"','').strip().split(",")
                ligand = linedata[ligandcol].upper()
                receptor = linedata[receptorcol].upper()
                if receptor not in lrdict.keys():
                    lrdict[receptor] = []
                lrdict[receptor].append(ligand)

    return lrdict

def loadCellTalkDB(lrfile, ligandcol=1, receptorcol=2):
    lrdict={}
    with open(lrfile,'r') as f:
        title = True
        for line in f.readlines():
            
            if title:
                title = False
            else:
                linedata = line.strip().split("\t")
                ligands = linedata[ligandcol].upper().strip().split('_')   
                receptors = linedata[receptorcol].upper().strip().split('_')                
                for receptor in receptors:
                    if receptor not in lrdict.keys():
                        lrdict[receptor] = []
                    for ligand in ligands:
                        lrdict[receptor].append(ligand)
    return lrdict


def splitrpacRoutes(rpacFile, pathways_folder, outputfile, do_write_to_file=False,
    ligand_list = ["LIGAND", "BUNDLE_LIGAND", "BUNDLE_OR_LIGAND","BUNDLE_AND_LIGAND", "BUNDLELIGANDAND", "BUNDLELIGANDOR"],
    receptor_list = ["RECEPTOR", "BUNDLE_RECEPTOR","BUNDLE_OR_RECEPTOR", "BUNDLE_AND_RECEPTOR"],
    tf_list = ["TRANSCRIPTION FACTOR", "BUNDLE_AND_TF", "BUNDLE_OR_TF","TRANSCRIPTIONFACTOR", "DIAMOND", "BUNDLETFCAND", "BUNDLETFCOR"],
):  
    
    line_sep='\n'
    route_col_sep='\t'
    header_row = ["ID","TYPE","LEN","PW1_NAME","PW1_ROUTE","PW1_ROUTE_GENES","P1_ligands","P1_receptors","TF"]
    all_route_scores = [header_row]

    if do_write_to_file:
        f_o = open(outputfile, "w")
        f_o.write(route_col_sep.join(header_row))

    all_routes = route.load_routes_from_file(rpacFile)
    all_routes_header_row = all_routes[0]
    graph_map = {}
    for route_row in all_routes[1:]:
        route_id = int(route_row[all_routes_header_row.index('ID')])
        
        # if(route_id==780):
        #     print("x")
        route_type = route_row[all_routes_header_row.index('TYPE')]
        pw1_name = route_row[all_routes_header_row.index('PW1_NAME')]
        if pw1_name not in graph_map:
            graph_map[pw1_name] = pathway_graph.load_graph_from_file(
                pathways_folder, pw1_name)
        G1 = graph_map[pw1_name]
        if G1 is None:
            continue
        pw1_route = route_row[all_routes_header_row.index('PW1_ROUTE')]
        # r1_e = score.get_o1_route_expectations(G1, pw1_route, route_type)
        P1_ligand = []        
        TF = []
        P1_receptor=[]
        
        this_route_scores = [
            'Route' + str(route_id),
            route_type,
            route_row[all_routes_header_row.index('LEN')],
            route_row[all_routes_header_row.index('PW1_NAME')],
            route_row[all_routes_header_row.index('PW1_ROUTE')],
            route_row[all_routes_header_row.index('PW1_ROUTE_GENES')]
            ]
        # if it is the P1 route
        # We have to sperate ligand and down stream
        if route_type =='P1':
            for node in pw1_route:
                if G1.nodes[node]['TYPE'].upper() in ligand_list:
                    P1_ligand.append(node)
                if G1.nodes[node]['TYPE'].upper() in receptor_list:
                    P1_receptor.append(node)
                if G1.nodes[node]['TYPE'].upper() in tf_list:
                   TF.append(node)
        
        if route_type =="P2":
            for node in pw1_route:
                if G1.nodes[node]['TYPE'].upper() in tf_list:
                    TF.append(node)

        this_route_scores.append(P1_ligand)
        this_route_scores.append(P1_receptor)
        this_route_scores.append(TF)

        if do_write_to_file:
            f_o.write(line_sep + route_col_sep.join([str(item)
                      for item in this_route_scores]))
        all_route_scores.append(this_route_scores)

    return all_route_scores,graph_map


def getnodesgenename(nids,G):
    genenames = []
    for nid in nids:    
        if not pathway_graph.is_bundle(G.nodes[nid]):
            genenames.append(G.nodes[nid]['NAME'])
        else:
            childrennodes  = G.graph['NODE_CHILDREN'][nid]
            genenames = genenames+getnodesgenename(childrennodes,G)
    return genenames





def getnodegenename(nid,G):
    genename=[]
    if not pathway_graph.is_bundle(G.nodes[nid]):
        genename.append(G.nodes[nid]['NAME'])
    else:
        childrennodes  = G.graph['NODE_CHILDREN'][nid]
        for childrenode in childrennodes:
            genename = genename+getnodegenename(childrenode,G)
    return genename


def getnodestype(nids,G):
    nodetypes = []
    for nid in nids:    
        if  pathway_graph.is_bundle_OR(G.nodes[nid]):
            nodetypes.append('OR_bundle')
        elif pathway_graph.is_bundle_AND(G.nodes[nid]):
            nodetypes.append('And_bundle')
        else:            
            nodetypes.append('gene')
    return nodetypes

def getnodetype(nid,G):
    
    
    if  pathway_graph.is_bundle_OR(G.nodes[nid]):
        return 'OR_bundle'
    elif pathway_graph.is_bundle_AND(G.nodes[nid]):
        return 'And_bundle'
    else: 
        return 'gene'


def generateCellCommunicationDB_V1(rpacRoutes,pathwaygraphmap,socketPairsDict,outputfile,do_write_to_file=False):
    # Go backwards
    all_routes_header_row = rpacRoutes[0]

    ResultsHeader = ['id','C2-R2-ID','C2-R2-PW','C2-R2-TYPE','C2-R2-ids','C2-R2-genes','C2-R2-RRecetors-ids','C2-R2-RRecetors-genes','C2-R2-RReceptors-type','C2-R2-RLigands-ids','C2-R2-RLigands-genes','C2-R2-RLigands-type','SocketPairs_Receptor','SocketPairs_Ligand','SocektPairs_TF','C1-R1-TFs-ids','C1-R1-TFs-genes','C1-R1-ids','C1-R1-genes','C1-R1-PW','C1-R1-ID','C1-R1-TYPE','combineType']

    '''
    id -> link id
    C2-R2 -> Route 2 ids 
    C2-R2-genes -> Route 2 genes
    C2-R2-RRecetors-ids -> Route 2 Receptors ids 
    C2-R2-RRecetors-genes -> Route 2 Receptors genes 
    C2-R2-RReceptors-type-> Route 2 Receptors type (gene/bundle) 
    C2-R2-RLigands-ids -> Route 2 Ligands ids 
    C2-R2-RLigands-genes -> Route 2 Ligands genes 
    C2-R2-RLigands-type-> Route 2 RLigands type (gene/bundle) 
    SocketPairs_Receptor -> Receptor in socket pairs
    SocketPairs_Ligand -> Ligand in socket pairs
    SocektPairs_TF -> TF in socket pairs
    C1-R1-TFs-ids -> Route 2 Ligands ids 
    C1-R1-TFs-geness -> Route 2 Ligands genes
    C1-R1 -> Route 1 ids 
    C1-R1-genes -> Route 1 genes
    

    '''


    # ResultsHeader = ['id','C1-R1','C1-R1-TF','C1-R1-TLigand','SocektPairs_TF','SocketPairs_Ligand','SocketPairs_Receptor','C2-R2-RLigand','C2-R2-RReceptor','C2-R2-R2']
    # socketPairs_header_row = ["ID","TF","Ligand(Target)","Receptor","Source"]
    ResultRows = [ResultsHeader]
    routepairs = {}
    index=1
    for R2_route_row in rpacRoutes[1:]:
        R2_route_id = R2_route_row[all_routes_header_row.index('ID')]
        print(R2_route_id)
        if (R2_route_id)=="Route780":
            print("x")
        R2_route_type = R2_route_row[all_routes_header_row.index('TYPE')]
        C2_R2_PW = R2_route_row[all_routes_header_row.index('PW1_NAME')]
        C2_R2_ids = R2_route_row[all_routes_header_row.index('PW1_ROUTE')]
        C2_R2_genes = R2_route_row[all_routes_header_row.index('PW1_ROUTE_GENES')]
        C2_R2_RReceptor_ids =  R2_route_row[all_routes_header_row.index('P1_receptors')]
        if R2_route_type == 'P1' and len(C2_R2_RReceptor_ids)>0:
            C2_R2_RLigand_ids = R2_route_row[all_routes_header_row.index('P1_ligands')]
            G2 = pathwaygraphmap[C2_R2_PW]
            C2_R2_RLigand_genes = getnodesgenename(C2_R2_RLigand_ids,G2)
            C2_R2_RLigand_type = getnodestype(C2_R2_RLigand_ids,G2)
            C2_R2_RReceptor_genes = getnodesgenename(C2_R2_RReceptor_ids,G2)
            C2_R2_RReceptor_type = getnodestype(C2_R2_RReceptor_ids,G2)
            if C2_R2_PW+R2_route_id not in routepairs.keys():
                routepairs[C2_R2_PW+R2_route_id] = []

            for R1_route_row in rpacRoutes[1:]:
                R1_route_id = R1_route_row[all_routes_header_row.index('ID')]
                R1_route_type = R1_route_row[all_routes_header_row.index('TYPE')]
                C1_R1_PW = R1_route_row[all_routes_header_row.index('PW1_NAME')]
                C1_R1_ids = R1_route_row[all_routes_header_row.index('PW1_ROUTE')]
                C1_R1_genes = R1_route_row[all_routes_header_row.index('PW1_ROUTE_GENES')]
                C1_R1_TFs_id =   R1_route_row[all_routes_header_row.index('TF')]
                G1 =  pathwaygraphmap[C1_R1_PW]
                C1_R1_TF_genes = getnodesgenename(C1_R1_TFs_id,G1)                    
                if R1_route_type == 'P2':
                    if len([i for i in C2_R2_RLigand_genes if i in C1_R1_genes])>0 and C1_R1_PW+R1_route_id not in routepairs[C2_R2_PW+R2_route_id]:
                        ResultRows.append(
                            [
                                str(index),
                                R2_route_id,
                                C2_R2_PW,
                                R2_route_type,
                            ','.join(C2_R2_ids),
                            ','.join(C2_R2_genes),
                            ','.join(C2_R2_RReceptor_ids),
                            ','.join(C2_R2_RReceptor_genes),
                            ','.join(C2_R2_RReceptor_type),
                            ','.join(C2_R2_RLigand_ids),                            
                            ','.join(C2_R2_RLigand_genes),
                            ','.join(C2_R2_RLigand_type),
                            ','.join(C2_R2_RReceptor_genes),
                            ','.join([i for i in C2_R2_RLigand_genes if i in C1_R1_genes]),
                            ','.join(C1_R1_TF_genes),
                            ','.join(C1_R1_TFs_id),
                            ','.join(C1_R1_TF_genes),
                            ','.join(C1_R1_ids),
                            ','.join(C1_R1_genes),
                            C1_R1_PW,
                            R1_route_id,R1_route_type,'2'])
                        index +=1
                        # print(index)
                        routepairs[C2_R2_PW+R2_route_id].append(C1_R1_PW+R1_route_id)
            for Receptor in C2_R2_RReceptor_genes:
                if Receptor in socketPairsDict.keys():
                    for ligand,tfs in socketPairsDict[Receptor].items(): # from socket pairs
                        if len(tfs)>0:
                            for tf in tfs: 
                                ResultRows.append(
                                    [str(index),
                                    R2_route_id,
                                    C2_R2_PW,
                                    R2_route_type,
                                    ','.join(C2_R2_ids),
                                    ','.join(C2_R2_genes),
                                    ','.join(C2_R2_RReceptor_ids),
                                    ','.join(C2_R2_RReceptor_genes),
                                    ','.join(C2_R2_RReceptor_type),
                                    ','.join(C2_R2_RLigand_ids),
                                    ','.join(C2_R2_RLigand_genes),
                                    ','.join(C2_R2_RLigand_type),
                                    Receptor,
                                    ligand,
                                    tf,
                                    'NA',
                                    'NA',
                                    'NA',
                                    'NA',
                                    'NA',
                                    'NA',
                                    'NA','4'])
                                index +=1
                                # print(index)
                        else:
                            ResultRows.append(
                                [str(index),
                                R2_route_id,
                                C2_R2_PW,
                                R2_route_type,
                                ','.join(C2_R2_ids),
                                ','.join(C2_R2_genes),
                                ','.join(C2_R2_RReceptor_ids),
                                ','.join(C2_R2_RReceptor_genes),
                                ','.join(C2_R2_RReceptor_type),
                                ','.join(C2_R2_RLigand_ids),
                                ','.join(C2_R2_RLigand_genes),
                                ','.join(C2_R2_RLigand_type),
                                Receptor,
                                ligand,
                                'NA',
                                'NA',
                                'NA',
                                'NA',
                                'NA',
                                'NA',
                                'NA',
                                'NA','5'])
                            index +=1
                            # print(index)


                            # for R1_route_row in rpacRoutes[1:]:
                            #     R1_route_id = R2_route_row[all_routes_header_row.index('ID')]
                            #     # R1_route_type = R2_route_row[all_routes_header_row.index('TYPE')]
                            #     C2_R2_PW = R2_route_row[all_routes_header_row.index('PW1_NAME')]
                            #     C1_R1_ids = R2_route_row[all_routes_header_row.index('PW1_ROUTE')]
                            #     C1_R1_genes = R2_route_row[all_routes_header_row.index('PW1_ROUTE_GENES')]
                            #     C1_R1_TFs_id =   R2_route_row[all_routes_header_row.index('TF')]
                            #     G1 =  pathwaygraphmap[R1_pw_name]
                            #     R1_TF_genes = getnodesgenename(C2_R2_RLigand,G1)
                            #     if tf in R1_TF_genes:

    if do_write_to_file:
        with open(outputfile,'w') as f:
            for row in ResultRows:
                f.write('\t'.join(row)+'\n')
    return ResultRows














# Here we difine the mode of the communication
# The communication should start from (C1-R1->C1-R1-TF->)C1-R1-TLigand -> C2-R2-Receptor -> C2-R2 



# def generateCellCommunicationDB_V2(rpacRoutes,pathwaygraphmap,socketPairsDict,outputfile,do_write_to_file=False):
#     # Go backwards
#     all_routes_header_row = rpacRoutes[0]

#     ResultsHeader = ['id','C2-R2-ID','C2-R2-PW','C2-R2-TYPE','C2-R2-ids','C2-R2-genes','C2-R2-RRecetors-ids','C2-R2-RRecetors-genes','C2-R2-RLigands-ids','C2-R2-RLigands-genes','SocketPairs_Receptor','SocketPairs_Ligand','SocektPairs_TF','C1-R1-TFs-ids','C1-R1-TFs-genes','C1-R1-ids','C1-R1-genes','C1-R1-PW','C1-R1-ID','C1-R1-TYPE','combineType']

#     '''
#     id -> link id
#     C2-R2 -> Route 2 ids 
#     C2-R2-genes -> Route 2 genes
#     C2-R2-RRecetors-ids -> Route 2 Receptors ids 
#     C2-R2-RRecetors-genes -> Route 2 Receptors genes 
#     C2-R2-RLigands-ids -> Route 2 Ligands ids 
#     C2-R2-RLigands-genes -> Route 2 Ligands genes 
#     SocketPairs_Receptor -> Receptor in socket pairs
#     SocketPairs_Ligand -> Ligand in socket pairs
#     SocektPairs_TF -> TF in socket pairs
#     C1-R1-TFs-ids -> Route 2 Ligands ids 
#     C1-R1-TFs-geness -> Route 2 Ligands genes
#     C1-R1 -> Route 1 ids 
#     C1-R1-genes -> Route 1 genes 
#     '''


#     # ResultsHeader = ['id','C1-R1','C1-R1-TF','C1-R1-TLigand','SocektPairs_TF','SocketPairs_Ligand','SocketPairs_Receptor','C2-R2-RLigand','C2-R2-RReceptor','C2-R2-R2']
#     # socketPairs_header_row = ["ID","TF","Ligand(Target)","Receptor","Source"]
#     ResultRows = [ResultsHeader]
#     routepairs = {}
#     index=1
#     for R2_route_row in rpacRoutes[1:]:
#         R2_route_id = R2_route_row[all_routes_header_row.index('ID')]
#         R2_route_type = R2_route_row[all_routes_header_row.index('TYPE')]
#         C2_R2_PW = R2_route_row[all_routes_header_row.index('PW1_NAME')]
#         C2_R2_ids = R2_route_row[all_routes_header_row.index('PW1_ROUTE')]
#         C2_R2_genes = R2_route_row[all_routes_header_row.index('PW1_ROUTE_GENES')]
#         if R2_route_type == 'P1':
#             C2_R2_RLigand_ids = R2_route_row[all_routes_header_row.index('P1_ligands')]
#             C2_R2_RReceptor_ids =  R2_route_row[all_routes_header_row.index('P1_receptors')]
#             G2 = pathwaygraphmap[C2_R2_PW]
#             C2_R2_RLigand_genes = getnodesgenename(C2_R2_RLigand_ids,G2)
#             C2_R2_RReceptor_genes = getnodesgenename(C2_R2_RReceptor_ids,G2)
#             if C2_R2_PW+R2_route_id not in routepairs.keys():
#                 routepairs[C2_R2_PW+R2_route_id] = []

#             for R1_route_row in rpacRoutes[1:]:
#                 R1_route_id = R1_route_row[all_routes_header_row.index('ID')]
#                 R1_route_type = R1_route_row[all_routes_header_row.index('TYPE')]
#                 C1_R1_PW = R1_route_row[all_routes_header_row.index('PW1_NAME')]
#                 C1_R1_ids = R1_route_row[all_routes_header_row.index('PW1_ROUTE')]
#                 C1_R1_genes = R1_route_row[all_routes_header_row.index('PW1_ROUTE_GENES')]
#                 C1_R1_TFs_id =   R1_route_row[all_routes_header_row.index('TF')]
#                 G1 =  pathwaygraphmap[C1_R1_PW]
#                 C1_R1_TF_genes = getnodesgenename(C1_R1_TFs_id,G1)                    
#                 for Receptor in C2_R2_RReceptor_genes:
#                     if Receptor in socketPairsDict.keys():
#                         for ligand,tfs in socketPairsDict[Receptor].items(): # from socket pairs
#                             for tf in tfs:
#                                 if tf in C1_R1_TF_genes:
#                                     ResultRows.append(
#                                         [str(index),
#                                         R2_route_id,
#                                         C2_R2_PW,
#                                         R2_route_type,
#                                         ','.join(C2_R2_ids),
#                                         ','.join(C2_R2_genes),
#                                         ','.join(C2_R2_RReceptor_ids),
#                                         ','.join(C2_R2_RReceptor_genes),
#                                         ','.join(C2_R2_RLigand_ids),
#                                         ','.join(C2_R2_RLigand_genes),
#                                         Receptor,
#                                         ligand,
#                                         tf,
#                                         ','.join(C1_R1_TFs_id),
#                                         ','.join(C1_R1_TF_genes),
#                                         ','.join(C1_R1_ids),
#                                         ','.join(C1_R1_genes),
#                                         C1_R1_PW,
#                                         R1_route_id,R1_route_type,'1'])
#                                     index +=1
#                                     print(index)
#                 if R1_route_type == 'P2':
#                     if len([i for i in C2_R2_RLigand_genes if i in C1_R1_genes])>0 and C1_R1_PW+R1_route_id not in routepairs[C2_R2_PW+R2_route_id]:
#                         ResultRows.append(
#                             [str(index),
#                             R2_route_id,
#                             C2_R2_PW,
#                             R2_route_type,
#                             ','.join(C2_R2_ids),
#                             ','.join(C2_R2_genes),
#                             ','.join(C2_R2_RReceptor_ids),
#                             ','.join(C2_R2_RReceptor_genes),
#                             ','.join(C2_R2_RLigand_ids),
#                             ','.join(C2_R2_RLigand_genes),
#                             ','.join(C2_R2_RReceptor_ids),
#                             ','.join(C2_R2_RReceptor_genes),
#                             ','.join([i for i in C2_R2_RLigand_genes if i in C1_R1_genes]),
#                             ','.join(C1_R1_TF_genes),
#                             ','.join(C1_R1_TFs_id),
#                             ','.join(C1_R1_TF_genes),
#                             ','.join(C1_R1_ids),
#                             ','.join(C1_R1_genes),
#                             C1_R1_PW,
#                             R1_route_id,R1_route_type,'2'])
#                         index +=1
#                         print(index)
#                         routepairs[C2_R2_PW+R2_route_id].append(C1_R1_PW+R1_route_id)

#                 if R1_route_type == 'P1':
#                     if C1_R1_PW+R1_route_id == C2_R2_PW+R2_route_id:
#                         ResultRows.append(
#                             [str(index),
#                             R2_route_id,
#                             C2_R2_PW,
#                             R2_route_type,
#                             ','.join(C2_R2_ids),
#                             ','.join(C2_R2_genes),
#                             ','.join(C2_R2_RReceptor_ids),
#                             ','.join(C2_R2_RReceptor_genes),
#                             ','.join(C2_R2_RLigand_ids),
#                             ','.join(C2_R2_RLigand_genes),
#                             'NA',
#                             'NA',
#                             'NA',
#                             'NA',
#                             'NA',
#                             'NA',
#                             'NA',
#                             C1_R1_PW,
#                             R1_route_id,R1_route_type,'3']
#                             )
#                         index +=1
#                         print(index)

#             for Receptor in C2_R2_RReceptor_genes:
#                     if Receptor in socketPairsDict.keys():
#                         for ligand,tfs in socketPairsDict[Receptor].items(): # from socket pairs
#                             if len(tfs)>0:
#                                 for tf in tfs: 
#                                     ResultRows.append(
#                                         [str(index),
#                                         R2_route_id,
#                                         C2_R2_PW,
#                                         R2_route_type,
#                                         ','.join(C2_R2_ids),
#                                         ','.join(C2_R2_genes),
#                                         ','.join(C2_R2_RReceptor_ids),
#                                         ','.join(C2_R2_RReceptor_genes),
#                                         ','.join(C2_R2_RLigand_ids),
#                                         ','.join(C2_R2_RLigand_genes),
#                                         Receptor,
#                                         ligand,
#                                         tf,
#                                         'NA',
#                                         'NA',
#                                         'NA',
#                                         'NA',
#                                         'NA',
#                                         'NA',
#                                         'NA','4'])
#                                     index +=1
#                                     print(index)
#                             else:
#                                 ResultRows.append(
#                                     [str(index),
#                                     R2_route_id,
#                                     C2_R2_PW,
#                                     R2_route_type,
#                                     ','.join(C2_R2_ids),
#                                     ','.join(C2_R2_genes),
#                                     ','.join(C2_R2_RReceptor_ids),
#                                     ','.join(C2_R2_RReceptor_genes),
#                                     ','.join(C2_R2_RLigand_ids),
#                                     ','.join(C2_R2_RLigand_genes),
#                                     Receptor,
#                                     ligand,
#                                     'NA',
#                                     'NA',
#                                     'NA',
#                                     'NA',
#                                     'NA',
#                                     'NA',
#                                     'NA',
#                                     'NA','5'])
#                                 index +=1
#                                 print(index)


#                                 # for R1_route_row in rpacRoutes[1:]:
#                                 #     R1_route_id = R2_route_row[all_routes_header_row.index('ID')]
#                                 #     # R1_route_type = R2_route_row[all_routes_header_row.index('TYPE')]
#                                 #     C2_R2_PW = R2_route_row[all_routes_header_row.index('PW1_NAME')]
#                                 #     C1_R1_ids = R2_route_row[all_routes_header_row.index('PW1_ROUTE')]
#                                 #     C1_R1_genes = R2_route_row[all_routes_header_row.index('PW1_ROUTE_GENES')]
#                                 #     C1_R1_TFs_id =   R2_route_row[all_routes_header_row.index('TF')]
#                                 #     G1 =  pathwaygraphmap[R1_pw_name]
#                                 #     R1_TF_genes = getnodesgenename(C2_R2_RLigand,G1)
#                                 #     if tf in R1_TF_genes:

#     if do_write_to_file:
#         with open(outputfile,'w') as f:
#             for row in ResultRows:
#                 f.write('\t'.join(row)+'\n')
#     return ResultRows








            # for socketpair_row in socketPairs[1:]:
            #     socketpair_TF = socketpair_row[socketPairs_header_row.index('TF')]
            #     socketpair_ligand = socketpair_row[socketPairs_header_row.index('Ligand(Target)')]
            #     socketpair_receptor = socketpair_row[socketPairs_header_row.index('Receptor')]
            #     socketpair_source = socketpair_row[socketPairs_header_row.index('Source')]
            #     if socketpair_receptor in R2_receptor_genes: 
                    # It meet the requirement of cell communication but we can go back to 



































# def paraseCellCommunicationDB(rpacRoutes,pathwaygraphmap,socketPairs,outputfile,do_write_to_file=False):
#     # allroutes = route.load_routes_from_file(rpacFile)
#     # The model we start from P2/P1 -> TF -> ligands -> receptors -> P1
#     #  header_row = ["ID","TYPE","LEN","PW1_NAME","PW1_ROUTE","PW1_ROUTE_GENES","P1_ligand","P1_receptors","P1_TF","P2_TF"]
#     #  => header_row = ["ID","PW1_TYPE","PW1_NAME","PW1_ROUTE","PW1_ROUTE_GENES","PW1_TF","PW1_target_ligand","PW2_ligand","PW2_recepotors","PW2_TYPE","PW2_NAME","PW2_DownRoute","PW2_ROUTE_GENES"]
#     all_routes_header_row = rpacRoutes[0]
#     header_row = ["ID",'R1_ID',"R1_TYPE","R1_NAME","R1_ROUTE","R1_ROUTE_GENES","R1_TF","R1_target_ligand",'socket_TF','socket_ligand','socket_receptor','R2_ID',"R2_ligand","R2_recepotors","R2_TYPE","R2_NAME","R2_Route","PW2_ROUTE_GENES","Source"]
#     socketPairs_header_row = ["ID","TF","Ligand(Target)","Receptor","Source"]

#     CommunicationModelroutes = [header_row]

#     for route_row in rpacRoutes[1:]:
#         # We start with each route if it P1 
#         R1_route_id = (route_row[all_routes_header_row.index('ID')])
#         R1_route_type = route_row[all_routes_header_row.index('TYPE')]
#         R1_pw_name = route_row[all_routes_header_row.index('PW1_NAME')]
#         R1_pw_route = route_row[all_routes_header_row.index('PW1_ROUTE')]
#         R1_pw_genes = route_row[all_routes_header_row.index('PW1_ROUTE_GENES')]
#         if len(route_row[all_routes_header_row.index('P1_TF')])>0:
#             R1_TF = route_row[all_routes_header_row.index('P1_TF')]
#         if len(route_row[all_routes_header_row.index('P2_TF')])>0:
#             R1_TF = route_row[all_routes_header_row.index('P2_TF')]

#         G1 = pathwaygraphmap[R1_pw_name]

#         R1_tf_genes = getnodesgenename(R1_TF,G1)
#         index = 1
#         for socketpair_row in socketPairs[1:]:
#             socketpair_TF = socketpair_row[all_routes_header_row.index('TF')]
#             socketpair_ligand = socketpair_row[all_routes_header_row.index('Ligand')]
#             socketpair_receptor = socketpair_row[all_routes_header_row.index('Receptor')]
#             socketpair_source = socketpair_row[all_routes_header_row.index('Source')]
#             if socketpair_TF in R1_tf_genes:
#                 for R2_route_row in rpacRoutes[1:]:
#                     # R2 should be P1
#                     R2_route_id = int(R2_route_row[all_routes_header_row.index('ID')])
#                     R2_route_type = R2_route_row[all_routes_header_row.index('TYPE')]
#                     if R2_route_type =='P1':
#                         R2_pw_name = R2_route_row[all_routes_header_row.index('PW1_NAME')]
#                         R2_pw_route = R2_route_row[all_routes_header_row.index('PW1_ROUTE')]
#                         R2_pw_ligand =  R2_route_row[all_routes_header_row.index('P1_ligand')]
#                         R2_pw_receptor =  R2_route_row[all_routes_header_row.index('P1_receptors')]
#                         R2_pw_genes = R2_route_row[all_routes_header_row.index('PW1_ROUTE_GENES')]
#                         G2 = pathwaygraphmap[R2_pw_name]
#                         R2_receptor_genes = getnodesgenename(R2_pw_receptor,G2)
#                         if socketpair_receptor in R2_receptor_genes:
#                             CommunicationModelroutes.append([str(index),R1_route_id,R1_route_type,R1_pw_name,R1_pw_route,R1_pw_genes,R1_TF,socketpair_TF,socketpair_ligand,socketpair_receptor,R2_route_id,R2_pw_ligand,R2_pw_receptor,R2_route_type,R2_pw_name,R2_pw_route,R2_pw_genes,socketpair_source])
#                             index+=1

#         if do_write_to_file:
#             with open(outputfile,'w') as f:
#                 for row in CommunicationModelroutes:
#                     f.write('\t'.join(row)+'\n')
#         return CommunicationModelroutes


        







        # for tfgene in R1_tf_genes:
        #     if tfgene in ChipSeqDB.keys():
        #         if len(targetlist)>0:
        #             targetlist=targetlist&ChipSeqDB[tfgene]
        #         else:
        #             targetlist=ChipSeqDB[tfgene]
        



            
        # if G1 is None:
        #     continue
        




# pathwaykgmlfolder='./BioKnowledge/KEGG/KGML/'
# pathwayjsonfolder='./BioKnowledge/KEGG/JSON/'
# Formpathwaycomponets(pathwaykgmlfolder,pathwayjsonfolder)







































