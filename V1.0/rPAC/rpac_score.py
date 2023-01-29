import sys, time
import json
from rPAC import rpac_pathway_graph as pathway_graph
from rPAC import rpac_util as util

NODE_CONTRIBUTION_MAX = 1
ZERO_NODE_PENALTY = 1

line_sep = "\n"
route_col_sep = "\t"

def set_node_contribution_max(node_contribution_max):
    global NODE_CONTRIBUTION_MAX
    NODE_CONTRIBUTION_MAX = node_contribution_max

def set_zero_node_penalty(zero_node_penalty):
    global ZERO_NODE_PENALTY
    ZERO_NODE_PENALTY = zero_node_penalty

def score_routes_from_gene_data_file(pathways_folder, all_routes, gene_data_file, gene_data_delimiter = "\t", max_sample_count=0, do_write_to_file=False, o_scores_file=""):
    header_row = all_routes[0]
    all_route_gene_set = set()
    genes_fields = ['PW1_ROUTE_GENES', 'ROUTE_GENES']
    for route in all_routes:
        for genes_field in genes_fields:
            if genes_field in header_row:
                all_route_gene_set.update(route[header_row.index(genes_field)])

    all_route_gene_set.add('GENE')
    gene_data = util.load_gene_data(
        gene_data_file, gene_data_delimiter, all_route_gene_set, max_sample_count)

    # print(gene_data)
    # print(len(gene_data['KOFAST_vs_WTFAST']))
    return score_routes(pathways_folder, all_routes, gene_data, do_write_to_file, o_scores_file)


def score_routes(pathways_folder, all_routes, gene_data, do_write_to_file=False, o_scores_file=""):
    all_routes_header_row = all_routes[0]
    # ['TYPE', 'LEN', 'PW1_NAME', 'PW1_ROUTE', 'PW1_ROUTE_NODE_EXP_LIST', 'PW1_ROUTE_GENES'
    # ,'PW1_SOCKET_GENES', 'NO_OF_XTALK_ROUTES', 'PW2_NAME', 'PW2_ROUTE', 'PW2_ROUTE_NODE_EXP_LIST', 'PW2_ROUTE_GENES', 'PW2_ROUTE_SOCKET_GENES']
    sample_list = gene_data['GENE']

    # header_row = all_routes_header_row + ['PW1_ROUTE_EXP_LIST'] + sample_list
    header_row = ['ID'] + sample_list
    if do_write_to_file:
        f_o = open(o_scores_file, "w")
        f_o.write(route_col_sep.join(header_row))

    all_route_scores = [header_row]
    graph_map = {}
    total_time_taken = 0
    start = time.process_time()
    route_count = 0
    prev_count = 0
    for route_row in all_routes[1:]:
        route_id = int(route_row[all_routes_header_row.index('ID')])
        # only run for certain routes
        # if route_id < 4652 or route_id > 4667:
        #     continue

        route_type = route_row[all_routes_header_row.index('TYPE')]
        pw1_name = route_row[all_routes_header_row.index('PW1_NAME')]

        if pw1_name not in graph_map:
            graph_map[pw1_name] = pathway_graph.load_graph_from_file(
                pathways_folder, pw1_name)
        G1 = graph_map[pw1_name]
        if G1 is None:
            continue

        pw1_route = route_row[all_routes_header_row.index('PW1_ROUTE')]

        # O1 routes
        if 'O2' not in route_type:
            r1_e = get_o1_route_expectations(
                G1, pw1_route, route_type)

            # this_route_scores = route_row.copy() + [r1_e]            
            this_route_scores = ['Route' + str(route_id)]
            # time_taken = 0
            for sample in sample_list:
                one_sample_gene_data = gene_data[sample]
                one_route_score = score_o1_route(
                    one_sample_gene_data, G1, pw1_route, r1_e)
                this_route_scores.append(one_route_score)
          
        if do_write_to_file:
            f_o.write(line_sep + route_col_sep.join([str(item)
                      for item in this_route_scores]))
        else:
            all_route_scores.append(this_route_scores)

        route_count += 1
        if route_count % 500 == 0:
            end = time.process_time()
            time_taken = end-start
            total_time_taken += time_taken
            print("routes", prev_count + 1, "to", route_count, "scored for",
                  len(sample_list), "samples in", round(total_time_taken, 2), "seconds")
            prev_count = route_count
            start = time.process_time()
    if route_count > prev_count:
        end = time.process_time()
        time_taken = end-start
        total_time_taken += time_taken
        print("routes", prev_count + 1, "to", route_count, "scored for",
                len(sample_list), "samples in", round(total_time_taken, 2), "seconds")
    # print("Total time taken for loading gene data and pathway ", total_time_taken)
    if do_write_to_file:
        f_o.close()
    return all_route_scores

def get_o1_route_expectations(G1, r1, route_type):
    r1_e = []
    r1_prev_nid = ''
    r1_nid_exp = 1
    r1_prev_nid_exp = 1
    if 'P1' in route_type:
        for r1_nid in r1[::-1]:
            if r1_prev_nid != '':
                r1_nid_exp = get_node_expectation(
                    G1, r1_nid, r1_prev_nid, r1_prev_nid_exp)
                r1_prev_nid_exp = r1_nid_exp
            r1_prev_nid = r1_nid
            r1_e.insert(0, r1_prev_nid_exp)

    elif 'P2' in route_type:
        for r1_nid in r1:
            if r1_prev_nid != '':
                # print(r1_prev_nid, r1_nid)
                r1_nid_exp = get_node_expectation(
                    G1, r1_prev_nid, r1_nid, r1_prev_nid_exp)
                r1_prev_nid_exp = r1_nid_exp
            r1_prev_nid = r1_nid
            r1_e.append(r1_prev_nid_exp)

    return r1_e


def get_raw_route_score(gene_data, G, r, r_e):
    zero_node_penalty_factor = ZERO_NODE_PENALTY  # ranges from 0 to 1
    raw_route_score = 0
    zero_nodes_count = 0
    r_len = len(r)
    route_len = r_len

    # scoring route
    for i in range(r_len):
        nid = r[i]
        node = G.nodes[nid]
        nid_exp = r_e[i]
        if pathway_graph.is_dummy_node(G.nodes[nid]):
            route_len -= 1
            continue
        node_values = get_node_values(gene_data, G, nid)
        extra_node_count = len(node_values) - 1
        route_len += extra_node_count
        for node_value in node_values:
            if node_value == 0:
                # zero_nodes_count += 1
                # route_len -= zero_nodes_count*(1 - zero_node_penalty_factor)                
                route_len -= (1 - zero_node_penalty_factor)
                # route_len -= 1
                continue
            else:
                node_contribution = get_node_contribution(node_value)
                node_score = node_contribution * nid_exp
                raw_route_score += node_score

    return (raw_route_score, route_len)


def score_o1_route(gene_data, G, r, r_e):
    route_score = 0
    (raw_route_score, route_len) = get_raw_route_score(
        gene_data, G, r, r_e)
    if route_len > 0 and raw_route_score != 0:
        route_score = round(raw_route_score/route_len, 2)
    if abs(route_score) > 1:
        route_score = route_score/abs(route_score)
    return route_score


def get_node_expectation(G, n1_id, n2_id, expectation):
    if G.edges[n1_id, n2_id]['STARTARROW'] == 'activate':
        return expectation
    elif G.edges[n1_id, n2_id]['STARTARROW'] == 'inhibit':
        return expectation * -1


def get_node_contribution(node_value):
    node_contribution_max = NODE_CONTRIBUTION_MAX
    node_contribution = node_value
    if abs(node_value) > node_contribution_max:
        node_contribution = node_contribution_max * node_value/abs(node_value)

    return round(node_contribution, 2)


def get_node_values(gene_data, G, nid):
    child_node_values = []
    node = G.nodes[nid]
    if pathway_graph.is_bundle_FG(node):
        child_node_ids = G.graph['NODE_CHILDREN'][nid]
        for child_node_id in child_node_ids:
            child_node_value = get_node_value(gene_data, G, child_node_id)
            if child_node_value != 0:
                child_node_values.append(child_node_value)
    else:
        node_value = get_node_value(gene_data, G, nid)
        child_node_values.append(node_value)
    return child_node_values


def get_node_value(gene_data, G, nid):
    node_value = 0
    node = G.nodes[nid]
    node_name = node['NAME']
    if not pathway_graph.is_bundle(node):
        if node_name in gene_data:
            node_val = gene_data[node_name]
            if is_number(node_val):
                node_value = float(node_val)
        # node_value = float(node['RNA'])
    elif pathway_graph.is_bundle_AND(node):
        node_value = evaluate_bundle_AND(gene_data, G, nid)
    elif pathway_graph.is_bundle_OR(node):
        node_value = evaluate_bundle_OR(gene_data, G, nid)
    elif pathway_graph.is_bundle_FG(node):
        node_value = evaluate_bundle_FG(gene_data, G, nid)
    else:
        node_value = evaluate_bundle_OR(gene_data, G, nid)
    return node_value


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def evaluate_bundle_OR(gene_data, G, nid):
    child_node_ids = G.graph['NODE_CHILDREN'][nid]

    child_node_values = []
    for child_node_id in child_node_ids:
        child_node_value = get_node_value(gene_data, G, child_node_id)
        if child_node_value != 0:
            child_node_values.append(child_node_value)

    if len(child_node_values) > 0:
        return max(child_node_values)
    return 0


def evaluate_bundle_AND(gene_data, G, nid):
    child_node_ids = G.graph['NODE_CHILDREN'][nid]
    child_node_values = []
    for child_node_id in child_node_ids:
        child_node_value = get_node_value(gene_data, G, child_node_id)
        if child_node_value != 0:
            child_node_values.append(child_node_value)

    if len(child_node_values) > 0:
        return min(child_node_values)
    return 0


def evaluate_bundle_FG(gene_data, G, nid):
    child_node_ids = G.graph['NODE_CHILDREN'][nid]
    child_node_values = []
    for child_node_id in child_node_ids:
        child_node_value = get_node_value(gene_data, G, child_node_id)
        if child_node_value != 0:
            child_node_values.append(child_node_value)
    if len(child_node_values) > 0:
        return max(child_node_values)
    return 0


def get_gene_data(G, nid, gene_data):
    gene_name = G.nodes[nid]['NAME']
    return gene_data[gene_name]


def load_gene_data_for_routes(gene_data_file, all_routes, max_sample_count):
    header_row = all_routes[0]
    all_route_gene_set = set()
    genes_fields = ['PW1_ROUTE_GENES', 'PW2_ROUTE_GENES', 'ROUTE_GENES']
    for route in all_routes:
        for genes_field in genes_fields:
            if genes_field in header_row:
                all_route_gene_set.update(route[header_row.index(genes_field)])

    all_route_gene_set.add('GENE')
    gene_data = util.load_gene_data(
        gene_data_file, all_route_gene_set, max_sample_count)

    return gene_data
