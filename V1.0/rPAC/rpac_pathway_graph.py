import sys, json
from numpy import ceil
import matplotlib.pyplot as plt
import networkx as nx
from rPAC import rpac_util as util
from os import path

def initGraph(json_file_path, gene_data={}):
    G = nx.DiGraph()
    global route_end_candidates
    json_obj = json.load(open(json_file_path))

    G.graph['NAME'] = json_obj['data']['NAME']
    node_attributes = {}
    node_children = {}
    for n in json_obj['elements']['nodes']:        
        this_node_attrs = {}
        for attr_key in n:
            if attr_key == 'data':
                for attr_key2 in n['data']:
                    this_node_attrs[attr_key2.upper()] = n['data'][attr_key2]
            else:
                this_node_attrs[attr_key.upper()] = n[attr_key]

        # if this node has parent, add this node as child of the parent
        if "PARENT" in this_node_attrs:
            parent_id = this_node_attrs["PARENT"]
            if parent_id.strip() == '':
                continue
            if parent_id not in node_children:
                node_children[parent_id] = []

            node_children[parent_id].append(this_node_attrs['ID'])

        this_node_id = this_node_attrs['ID']
        this_node_name = this_node_attrs['NAME']
        this_node_type = this_node_attrs['TYPE'].upper()
        #also load gene data if it is provided
        this_node_attrs['RNA'] = 0
        if this_node_name in gene_data:
            this_gene_data = gene_data[this_node_name]
            if util.is_number(this_gene_data):
                this_node_attrs['RNA'] = this_gene_data

        G.add_node(this_node_id)
        node_attributes[this_node_id] = this_node_attrs

    nx.set_node_attributes(G, node_attributes)

    for e in json_obj['elements']['edges']:
        this_edge_source = e['data']['source']
        this_edge_target = e['data']['target']
        this_edge_type = e['data']['StartArrow']
        if this_edge_type not in ['activate', 'inhibit']:
            continue
        this_edge_attrs = {}
        for edge_attr_key in e:
            if edge_attr_key == 'data':
                for attr_key in e['data']:
                    this_edge_attrs[attr_key.upper()] = e['data'][attr_key]
            else:
                this_edge_attrs[edge_attr_key.upper()] = e[edge_attr_key]

        G.add_edges_from(
            [(this_edge_source, this_edge_target, this_edge_attrs)])

    G.graph['NODE_CHILDREN'] = node_children
    # print(G.graph['NODE_CHILDREN'])
    # sys.exit()
    return G

def load_gene_data_to_graph(G, gene_data):
    for nid in G.nodes:
        this_node_name = G.nodes[nid]['NAME']
        # load gene data
        G.nodes[nid]['RNA'] = 0
        if this_node_name in gene_data:
            this_gene_data = gene_data[this_node_name]
            if util.is_number(this_gene_data):
                G.nodes[nid]['RNA'] = this_gene_data
    return G

def find_nodes_for_genes(G, genes = []):
    result_nodes = {}
    nodes = G.nodes
    leaf_nodes = list(v for v, d in G.out_degree() if d == 0)
    for nid in nodes:
        # exclude leaf nodes
        if nid in leaf_nodes:
            continue
        this_node_genes = get_node_genes(G, nid)
        for gene in this_node_genes:
            if len(gene.strip()) == 0:
                continue
            if gene in genes or len(genes) == 0:
                if gene not in result_nodes:
                    result_nodes[gene] = []
                if nid not in result_nodes[gene]:
                    result_nodes[gene].append(nid)
    return result_nodes

def find_gene_nodes_for_pathway(json_file_path):
    G = initGraph(json_file_path)
    nodes = find_nodes_for_genes(G)
    return nodes

def find_unique_nodes_for_genes(G, genes):
    result_nodes = {}
    nodes = G.nodes
    for nid in nodes:
        this_node_genes = get_node_genes(G, nid)
        for gene in this_node_genes:
            if gene in genes:
                if nid not in result_nodes:
                    # result_nodes.append(nid)
                    result_nodes[nid] = []
                result_nodes[nid].append(gene)

    return result_nodes


def is_pathway_node(node):
    if node['TYPE'] == 'PATHWAY':
        return True

    return False


def is_dummy_node(n):
    singleton_node_types = ['LIGAND', 'RECEPTOR',
                            'TRANSCRIPTION FACTOR', 'GENE']

    node_type = n['TYPE'].upper()
    if node_type in singleton_node_types or is_bundle(n):
        return False

    return True


def is_bundle(n):
    node_type = n['TYPE'].upper()
    if node_type.find('BUNDLE') >= 0:
        return True
    return False


def is_bundle_FG(n):
    node_type = n['TYPE'].upper()
    if node_type == 'BUNDLE' or node_type == 'BUNDLE_FG':
        return True
    return False


def is_bundle_OR(n):
    node_type = n['TYPE'].upper()
    if node_type.find('OR') >= 0:
        return True
    return False


def is_bundle_AND(n):
    node_type = n['TYPE'].upper()
    if node_type.find('AND') >= 0:
        return True
    return False



def get_node_genes(G, nid):
    genes = []
    node = G.nodes[nid]    
    if nid not in G.graph['NODE_CHILDREN']:
        genes = [node['NAME']]
    else:
        for child_nid in G.graph['NODE_CHILDREN'][nid]:
            genes += get_node_genes(G, child_nid)

    return genes

def get_node_nodes(G, nid):    
    nodes = [nid]
    if nid in G.graph['NODE_CHILDREN']:
        for child_nid in G.graph['NODE_CHILDREN'][nid]:
            nodes += get_node_nodes(G, child_nid)

    return nodes

def visualize_routes(all_routes, i_pathway_folder, o_output_folder):
    nrows = 1
    ncols = 2
    route_no = 1
    for pw1 in all_routes:

        print(pw1)
        i_pathway_file_path = i_pathway_folder + pw1 + ".json"
        G1 = initGraph(i_pathway_file_path)
        routes_in_pathway = all_routes[pw1]

        for route_type in routes_in_pathway:
            routes = routes_in_pathway[route_type]
            for route in routes:
                pw1_route_nodes = route['PW1_ROUTE']
                socket_genes = route['SOCKET_GENES']
                pw2 = route['PW2_NAME']
                pw2_file_path = i_pathway_folder + "KEGG - " + pw2 + ".json"
                G2 = initGraph(pw2_file_path)
                pw2_routes_dict = route['PW2_ROUTES']
                # pw2_routes_count = len(pw2_routes)
                # total_routes_count = pw2_routes_count + 1
                # nrows = round(total_routes_count ** 0.5)

                for nid in pw2_routes_dict:
                    # print(routes_in_pathway)
                    pw2_this_node_routes = pw2_routes_dict[nid]['ROUTES']
                    this_node_socket_genes = pw2_routes_dict[nid]['ROUTE_SOCKET_GENES']
                    for pw2_route in pw2_this_node_routes:

                        figure_no = 1
                        fig = plt.figure(figsize=(12, 6))
                        ax = plt.gca()
                        fig.text(0.5, -0.13, "Crosstalk genes " + str(this_node_socket_genes), fontsize=15, transform=ax.transAxes, verticalalignment='bottom',
                                 horizontalalignment='center')

                        plt.subplot(nrows, ncols, figure_no).set_title(pw1)
                        figure_no += 1

                        draw_pathway(G1, pw1_route_nodes)
                        plt.subplot(nrows, ncols, figure_no).set_title(pw2)
                        figure_no += 1
                        draw_pathway(G2, pw2_route)
                        plt.subplots_adjust(left=0.05,
                                            bottom=0.05,
                                            right=0.95,
                                            top=0.95,
                                            wspace=0.05,
                                            hspace=0.05)
                        o2_route_name = pw1 + "_" + route_type + \
                            "_X_" + pw2 + "_" + str(route_no)
                        o2_routes_file_path = o_output_folder + o2_route_name

                        plt.savefig(o2_routes_file_path)
                        plt.clf()
                        plt.cla()
                        plt.close()

                        # sys.exit()
                        route_no += 1
                    pass

    pass


def draw_pathway(G, node_list=[], edge_list=[]):
    full_node_list = [node for node in node_list]
    pos = {}
    labels = {}
    nodes = G.nodes
    # prev_nid = ""
    node_colors = []
    all_labels = {}
    # print(G.graph['NAME'],  G.graph['NODE_CHILDREN'] )
    for nid in nodes:
        node = G.nodes[nid]
        pos[nid] = [node['POSITION']["x"], (-1) * node['POSITION']["y"]]
        all_labels[nid] = nodes[nid]['NAME']
        if nid in node_list:
            # node_colors.append('yellow')
            labels[nid] = nodes[nid]['NAME']
            child_nodes = get_node_nodes(G, nid)
            for child_nid in child_nodes:
                labels[child_nid] = nodes[child_nid]['NAME']
                full_node_list.append(child_nid)

        #     node_colors.append('#1f78b4')
        # if prev_nid != '':
        #     edges.append((prev_nid, nid))

        # prev_nid = nid

    if len(node_list) == 0:
        nodes = list(G.nodes)

    if len(edge_list) == 0:
        edges = list(G.edges())

    # nx.draw_networkx(G, labels=labels, pos=pos, nodelist=nodes,
    #                  edgelist=edges, node_color=node_colors, alpha = 0.2, font_color='grey', node_size = 100, font_size = 8)

    nx.draw_networkx(G, pos, labels=all_labels, alpha=0.2,
                     font_color='black', node_size=100, font_size=8)
    nx.draw_networkx_nodes(G, pos, nodelist=full_node_list, node_color='yellow')
    # print(labels)
    nx.draw_networkx_labels(G, pos=pos, labels=labels)


def visualize_routes_on_tiles(all_routes, i_pathway_folder, o_output_folder):
    nrows = 1
    ncols = 2
    route_no = 0
    for pw1 in all_routes:

        print(pw1)
        i_pathway_file_path = i_pathway_folder + pw1 + ".json"
        G1 = initGraph(i_pathway_file_path)
        routes_in_pathway = all_routes[pw1]

        

        for route_type in routes_in_pathway:
            routes = routes_in_pathway[route_type]
            for route in routes:
                
                pw1_route_nodes = route['PW1_ROUTE']
                socket_genes = route['SOCKET_GENES']
                pw2 = route['PW2_NAME']
                pw2_file_path = i_pathway_folder + "KEGG - " + pw2 + ".json"
                G2 = initGraph(pw2_file_path)
                pw2_routes_dict = route['PW2_ROUTES']
                # pw2_routes_count = len(pw2_routes)
                # total_routes_count = pw2_routes_count + 1
                # nrows = round(total_routes_count ** 0.5)

                 
                no_of_xtalk_routes = route['NO_OF_XTALK_ROUTES']
                # skip for routes without any cross talk
                if no_of_xtalk_routes == 0 or no_of_xtalk_routes > 16:
                    continue

                total_no_of_plots = no_of_xtalk_routes + 1
                nrows = 1
                ncols = 2
                
                sqrt = ceil(total_no_of_plots ** 0.5)
                no_of_plots_possible = 1
                while nrows * ncols < total_no_of_plots:
                    
                    if nrows == ncols:
                        ncols += 1
                    elif nrows < ncols:
                        nrows += 1

                    pass
                # row_col_map = {'2':(1,2), '3':(2,2), '4': (2,2), '5': (2,3), '6': (2,3), '7': (3,3), '8': (3,3), '9': (3,3)}

                # if total_no_of_plots == 2:
                #     nrows = 1
                #     ncols = 2
                
                # elif total_no_of_plots == 6:
                #     nrows = 2
                #     ncols = 3
                # else:
                #     nrows = ceil(total_no_of_plots ** 0.5)
                #     ncols = nrows

                fig = plt.figure(figsize=(12, 6))
                ax = plt.gca()
                
                figure_no = 1                
                plt.subplot(nrows, ncols, figure_no).set_title(pw1)
                figure_no += 1
                draw_pathway(G1, pw1_route_nodes)

                for nid in pw2_routes_dict:
                    # print(routes_in_pathway)
                    pw2_this_node_routes = pw2_routes_dict[nid]['ROUTES']
                    this_node_socket_genes = set(pw2_routes_dict[nid]['ROUTE_SOCKET_GENES'])
                    for pw2_route in pw2_this_node_routes:

                        route_no += 1

                        # figure_no = 1
                        # fig = plt.figure(figsize=(12, 6))
                        # ax = plt.gca()
                        # fig.text(0.5, -0.13, "Crosstalk genes " + str(this_node_socket_genes), fontsize=15, transform=ax.transAxes, verticalalignment='bottom',
                        #          horizontalalignment='center')
                        plt.suptitle(pw1 + " --X-- " + pw2) 
                        plt.subplot(nrows, ncols, figure_no).set_title(this_node_socket_genes)
                        figure_no += 1
                        draw_pathway(G2, pw2_route)
                        plt.subplots_adjust(left=0.05,
                                            bottom=0.05,
                                            right=0.95,
                                            top=0.9,
                                            wspace=0.05,
                                            hspace=0.2)
                       

                        # sys.exit()
                        
                o2_route_name = pw1 + "_" + route_type + \
                            "_X_" + pw2 + "_" + str(route_no)
                o2_routes_file_path = o_output_folder + o2_route_name

                plt.savefig(o2_routes_file_path)
                plt.clf()
                plt.cla()
                plt.close()
                # sys.exit()

    pass


def visualize_route(pathways_folder, route, header_row, o_image_file):
    # header_row = ['ID', 'TYPE', 'LEN', 'PW1_NAME', 'PW1_ROUTE',  'PW1_ROUTE_GENES',
    #               'PW1_SOCKET_GENES', 'NO_OF_XTALK_ROUTES', 'PW2_NAME', 'PW2_ROUTE', 'PW2_ROUTE_GENES', 'PW2_ROUTE_SOCKET_GENES']
    field_indices = {}
    for i in range(len(header_row)):
        field_indices[header_row[i]] = i
    
    route_id = route[field_indices['ID']]
    route_type = route[field_indices['TYPE']]
    pw1_name = route[field_indices['PW1_NAME']]
    pw1_route = route[field_indices['PW1_ROUTE']]
    
    pw2_name = ""
    if len(route) > field_indices['PW2_ROUTE']:    
        pw2_name = route[field_indices['PW2_NAME']]
        pw2_route = route[field_indices['PW2_ROUTE']]

    nrows = 1
    ncols = 1
    if 'O2' in route_type:
        ncols = 2
    
    G1 = load_graph_from_file(pathways_folder, pw1_name)    

    fig = plt.figure(figsize=(12, 6))
    ax = plt.gca()
    
    figure_no = 1                
    plt.subplot(nrows, ncols, figure_no).set_title(pw1_name, fontsize=25)
    figure_no += 1
    draw_pathway(G1, pw1_route)

    if ncols > 1:
        plt.subplot(nrows, ncols, figure_no).set_title(pw2_name, fontsize=25)
        G2 = load_graph_from_file(pathways_folder, pw2_name)
        draw_pathway(G2, pw2_route)
        plt.subplots_adjust(left=0.05,
                            bottom=0.05,
                            right=0.95,
                            top=0.9,
                            wspace=0.05,
                            hspace=0.2)
            

            # sys.exit()        

    plt.savefig(o_image_file)
    plt.clf()
    plt.cla()
    plt.close()
    # sys.exit()

    pass

def find_gene_nodes_for_pathways(i_pathway_folder):
    gene_nodes = {}
    pathway_files = util.get_all_pathway_files(i_pathway_folder)
    for pathway in pathway_files:
        pathway_file_path = i_pathway_folder + pathway + ".json"
        this_pathway_gene_nodes = find_gene_nodes_for_pathway(pathway_file_path)
        for gene in this_pathway_gene_nodes:
            nodes = this_pathway_gene_nodes[gene]
            if gene not in gene_nodes:
                gene_nodes[gene] = {}
            
            gene_nodes[gene][pathway] = nodes

    return gene_nodes

def get_route_end_candidates(G, route_end_candidate_types):
    # print(route_end_candidate_types)
    leaf_nodes = list(v for v, d in G.out_degree() if d == 0)
    
    end_types = ['SOURCE', 'TARGET']
    # initializing route_end_candidate_nodes dictionary
    route_end_candidate_nodes = {}
    end_node_type_list = []
    for route_type in route_end_candidate_types:
        route_end_candidate_nodes[route_type] = {}
        # end_type = 'SOURCE' 
        for end_type in end_types:
            if end_type in route_end_candidate_types[route_type]:
                route_end_candidate_nodes[route_type][end_type] = []
                end_node_type_list += route_end_candidate_types[route_type][end_type]
        
    if len(end_node_type_list) == 0:
        return route_end_candidate_nodes

    nodes = G.nodes
    for nid in nodes:
        this_node_id = nid
        this_node_type = nodes[nid]['TYPE'].upper()
        # if this_node_id == 'n369':
        #     print(this_node_type)
        
        for route_type in route_end_candidate_types:
            for end_type in end_types:
                if end_type in route_end_candidate_types[route_type]:
                    if this_node_type in route_end_candidate_types[route_type][end_type]:    
                        route_end_candidate_nodes[route_type][end_type].append(
                            this_node_id)
        
    for route_type in route_end_candidate_types:
        # post processing for source candidate nodes
        # remove source node candidate if it is a leaf node
        # in P1 route types, remove source node candidate if the upstream is also a source candidate
        end_type = 'SOURCE'
        if end_type not in route_end_candidate_nodes[route_type]:
            continue
        sources = [
            source for source in route_end_candidate_nodes[route_type][end_type]]
        for source in sources:
            # remove source node candidate if it is a leaf node
            if source in leaf_nodes:
                route_end_candidate_nodes[route_type][end_type].remove(
                    source)
        
        # in P1 route types, remove source node candidate if the upstream is also a source candidate
        if route_type == 'P1':
            edges = G.edges(sources)
            for (source, target) in edges:
                if target in route_end_candidate_nodes[route_type][end_type]:
                    # remove source node candiate if the upstream is also a source candidate
                    if source in route_end_candidate_nodes[route_type][end_type]:
                        route_end_candidate_nodes[route_type][end_type].remove(
                            target)

        end_type = 'TARGET'
        if end_type not in route_end_candidate_nodes[route_type]:
            continue
        if len(route_end_candidate_nodes[route_type][end_type]) == 0:
            leaf_allowed = False
            if 'LEAF_ALLOWED' in route_end_candidate_types[route_type]:
                leaf_allowed = route_end_candidate_types[route_type]['LEAF_ALLOWED']
            
            if leaf_allowed:
                route_end_candidate_nodes[route_type][end_type] = leaf_nodes
            else:
                continue

    return route_end_candidate_nodes

def load_graph_from_file(pathway_folder, pw_name):
    pw_file_name = pw_name + '.json'
    if not "KEGG - " in pw_file_name:
        pw_file_name = 'KEGG - ' + pw_file_name
    pathway_file_path = pathway_folder + pw_file_name
    if not path.exists(pathway_file_path):
        return None
    return initGraph(pathway_file_path)