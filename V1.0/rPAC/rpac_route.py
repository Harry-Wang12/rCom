
import networkx as nx
import ast
from rPAC import rpac_pathway_graph as pathway_graph
from rPAC import rpac_util as util

route_min_length = {'P1': 3, 'P2': 2}
route_col_sep = "\t"
line_sep = "\n"
def find_routes_on_pathway(G, route_end_candidate_types):
    route_end_candidate_nodes = pathway_graph.get_route_end_candidates(
        G, route_end_candidate_types)

    pathway = G.graph['NAME']
    pw1_name = pathway.replace("KEGG - ", "")

    all_o1_route_candidates = {}
    all_o1_routes = []

    for route_type in route_end_candidate_nodes:
        if 'SOURCE' not in route_end_candidate_nodes[route_type]:
            continue
        all_o1_route_candidates[route_type] = []
        this_route_type_end_candidates = route_end_candidate_nodes[route_type]
        source_candidate_nodes = this_route_type_end_candidates['SOURCE']
        target_candidate_nodes = this_route_type_end_candidates['TARGET']

        for source in source_candidate_nodes:
            routes = list(nx.all_simple_paths(
                G, source=source, target=target_candidate_nodes))

            for route in routes:
                this_route = []
                route_len = 0
                i = 0
                pw1_route_genes = []
                for nid in route:
                    if not pathway_graph.is_dummy_node(G.nodes[nid]):
                        route_len += 1
                    i += 1
                    pw1_route_genes += pathway_graph.get_node_genes(
                        G, nid)
                    if nid in target_candidate_nodes:
                        break

                if i == len(route) and route_len >= route_min_length[route_type]:
                    this_route = [route_type, route_len,
                                  pw1_name, route, pw1_route_genes]
                    all_o1_routes.append(this_route)

    return all_o1_routes


def find_routes_on_pathways(i_pathway_folder, route_end_candidate_types, do_write_to_file=False, o_o1_routes_file=''):
    pathway_files = util.get_all_pathway_files(i_pathway_folder)
    all_pathway_o1_routes = []
    header_row = ['ID', 'TYPE', 'LEN', 'PW1_NAME', 'PW1_ROUTE', 'PW1_ROUTE_GENES']
    all_pathway_o1_routes.append(header_row)
    route_id = 0
    for pathway in pathway_files:
        pathway_file_path = i_pathway_folder + pathway + ".json"
        # find routes for this particular pathway
        G = pathway_graph.initGraph(pathway_file_path)
        this_pathway_o1_routes = find_routes_on_pathway(
            G, route_end_candidate_types)
        for one_o1_route in this_pathway_o1_routes:
            route_id += 1
            all_pathway_o1_routes.append([route_id] + one_o1_route)
        
    if do_write_to_file:
        write_routes_to_file(all_pathway_o1_routes, o_o1_routes_file)

    return all_pathway_o1_routes

def write_routes_to_file(all_routes, o_routes_file_path):
    f_o = open(o_routes_file_path, "w")
    first_line = True
    new_line = ""
    for route in all_routes:
        route_row = [str(r) for r in route]
        line_to_write = new_line + route_col_sep.join(route_row)
        f_o.write(line_to_write)

        if first_line:
            new_line = line_sep
            first_line = False

    f_o.close()

def load_routes_from_file(route_file_path):
    f = open(route_file_path)
    route_list = []
    first_row = True
    header_row = []
    indices_to_eval = []
    for line in f:
        line_list = line.strip().split(route_col_sep)
        if first_row:
            header_row = line_list
            first_row = False
            one_route = line_list
            # The values in these indices need to be evaluated so proper data type is retained (eg. list and set)
            header_fields_to_eval = ['PW1_ROUTE', 'PW1_ROUTE_GENES']
            indices_to_eval = []
            for header_field in header_fields_to_eval:
                if header_field in header_row:
                    indices_to_eval.append(header_row.index(header_field))
        else:
            one_route = line_list            
            for i in indices_to_eval:
                val = line_list[i]
                if val != '':
                    one_route[i] = ast.literal_eval(line_list[i])

        route_list.append(one_route)
    
    f.close()
    return route_list
    