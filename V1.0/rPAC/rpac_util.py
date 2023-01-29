import os, json, sys, ast
route_col_sep = "\t"
# col_sep = "\t"
line_sep = "\n"
gene_data_delimiter = "\t"
# def set_delimiter(delimiter = '\t'):
#     global gene_data_delimiter
#     gene_data_delimiter = delimiter
def set_gene_data_delimiter(delimiter = '\t'):
    global gene_data_delimiter
    gene_data_delimiter = delimiter
# def write_routes_to_file(all_routes, o_routes_file_path):
#     f_o = open(o_routes_file_path, "w")
#     first_line = True
#     new_line = ""
#     for route in all_routes:
        
#         route_row = [str(r) for r in route]
#         line_to_write = new_line + route_col_sep.join(route_row)
#         f_o.write(line_to_write)

#         if first_line:
#             new_line = line_sep
#             first_line = False
#         pass

#     f_o.close()
#     pass

# def write_data_matrix_to_file(data_matrix, o_file_path):
#     f_o = open(o_file_path, "w")
#     first_line = True
#     new_line = ""
#     for row in data_matrix:
#         route_row = [str(col) for col in row]
#         line_to_write = new_line + route_col_sep.join(route_row)
#         f_o.write(line_to_write)

#         if first_line:
#             new_line = "\n"
#             first_line = False
#         pass

#     f_o.close()
#     pass

# def write_routes_and_scores_to_file(all_routes, all_route_scores, o2_route_score_file_path):
#     f_o = open(o2_route_score_file_path, "w")
#     first_line = True
#     route_count = len(all_routes)
    
#     for i in range(route_count):
#         route_dict = all_routes[i]
#         route_scores = all_route_scores[i+1]
#     # for route_dict in all_routes:
#         if first_line:

#             header_row_list = list(route_dict.keys())

#             for sample_name in all_route_scores[0]:
#                 header_row_list.append(sample_name)
#             header_line = route_col_sep.join(header_row_list)
#             f_o.write(header_line)
#             first_line = False
#         route_row_list = []
#         for key, val in route_dict.items():
#             route_row_list.append(str(val))
        
#         for route_score in route_scores:
#             route_row_list.append(str(route_score))

#         line_to_write = line_sep
#         line_to_write += route_col_sep.join(route_row_list)
#         f_o.write(line_to_write)
#         pass

#     f_o.close()
#     pass


def get_all_pathway_files(i_pathway_folder):
    pathway_files = []
    for file in os.listdir(i_pathway_folder):
        if file.endswith(".json"):
            # i_pathway_file_path = i_pathway_folder + file
            pathway_files.append(file.replace(".json", ""))
    return pathway_files


def load_conf(conf_file_path):
    json_obj = json.load(open(conf_file_path))
    return json_obj


def load_gene_data(i_data_file, gene_data_delimiter = "\t",  genes = [], max_sample_count = 0):
    gene_data = {}
    f = open(i_data_file)
    sample_list = []
    first_line = True
    for line in f:
        gene = line[:line.find(gene_data_delimiter)].upper()
        if first_line:
            gene = 'GENE'
            first_line = False
        if gene in genes or len(genes) == 0:
            line_list = line.strip().split(gene_data_delimiter)[1:]
            if max_sample_count > 0 and max_sample_count < len(line_list):
                line_list = line_list[:max_sample_count]
            if gene == 'GENE':
                sample_list = line_list
                gene_data[gene] = sample_list
            else:
                for i in range(len(sample_list)):
                    sample = sample_list[i]
                    if sample not in gene_data:
                        gene_data[sample] = {}
                    gene_data[sample][gene] = line_list[i]
                pass
    
    return gene_data

def load_gene_data_for_genes(gene_data, genes):
    new_gene_data = {'GENE': gene_data['GENE']}
    
    for gene in gene_data:
        if gene in genes:
            new_gene_data[gene] = gene_data[gene]

    return new_gene_data

def load_gene_data_for_sample(gene_data, sample_index ):
    new_gene_data = {}
    # sample_list = gene_data['GENE']
    # sample_index = sample_list.index(sample_name)
    for gene in gene_data:
        new_gene_data[gene] = gene_data[gene][sample_index]

    return new_gene_data

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# def load_routes_from_file(route_file_path):
#     f = open(route_file_path)
#     route_list = []
#     first_row = True
#     header_row = []
#     indices_to_eval = []
#     for line in f:
#         line_list = line.strip().split(route_col_sep)
#         if first_row:
#             header_row = line_list
#             first_row = False
#             one_route = line_list
#             # The values in these indices need to be evaluated so proper data type is retained (eg. list and set)
#             header_fields_to_eval = ['PW1_ROUTE', 'PW1_ROUTE_GENES']
#             indices_to_eval = []
#             for header_field in header_fields_to_eval:
#                 if header_field in header_row:
#                     indices_to_eval.append(header_row.index(header_field))
#         else:
#             one_route = line_list            
#             for i in indices_to_eval:
#                 val = line_list[i]
#                 if val != '':
#                     one_route[i] = ast.literal_eval(line_list[i])

#         route_list.append(one_route)
    
#     f.close()
#     return route_list
        