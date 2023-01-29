import rpac_route as route, rpac_score as score
import random, sys, time

# o_image_folder = "G:/My Drive/UCONN/Research/20210523SUMM21_research/cross_talk/analysis/rpax_analysis_20210825/routes/"

FIND_ROUTES = True
LOAD_ROUTES_FROM_FILE = False
SCORE_ROUTES = True
DRAW_ROUTES = False
suffix = time.strftime("%Y%m%d%H%M%S") + "_" + str(random.randint(1, 100))
max_sample_count = 0
node_contribution_max = 1
zero_node_penalty = 1
gene_data_delimiter = '\t'
suffix = 'NULL'
##### DO NOT CHANGE ANYTHING BEYOND THIS LINE #####
help_msg = '''
command line usage: 
python rpac.py
-I_PATHWAYS_FOLDER <folder with all pathways JSON files> (required)
-FIND_ROUTES <True|False> (default:False)
-O_ROUTES_FILE <output file to write routes> (required if finding routes)
-LOAD_ROUTES_FROM_FILE <True|False> (default:False)
-I_ROUTES_FILE <input file with routes> (required if loading routes from file)
-I_GENE_DATA_FILE <gene rna data file to use for scoring> (required for scoring)
-MAX_SAMPLE_COUNT <number of samples to use for scoring, default 0 for all samples>
-SCORE_ROUTE <True|False> (default:False)
-O_SCORES_FILE <output file to write scores> (required for scoring)
-SUFFIX <suffix to use for output file name, NULL for no suffix> (default: current datetime and random number) 
-ZERO_NODE_PENALTY <penalty weight for nodes in route with 0 values.> (default:0, ranges from 0 to 1)
-GENE_DATA_DELIMITER <comma|tab> (default:tab)
'''
# def parse_sys_argv(help_msg=help_msg):
#     sys_args = sys.argv
#     args_len  = len(sys_args)
#     if args_len < 2:
#         print(help_msg)
#     args = {} 
#     for i in range(len(sys_args)):
#         arg = sys_args[i]
#         if '-' in arg and len(sys_args) > (i+1):
#             args[arg] = sys_args[i+1]            
#     return args
# args = parse_sys_argv()
# print(args)
# if '-I_PATHWAYS_FOLDER' in args:
#     i_pathways_folder = args['-I_PATHWAYS_FOLDER']
# if '-FIND_ROUTES' in args:
#     if args['-FIND_ROUTES'].upper() == 'TRUE':
#         FIND_ROUTES = True
# if '-SUFFIX' in args:
#     suffix = args['-SUFFIX'].strip()
# if '-O_ROUTES_FILE' in args:
#     o_routes_file = args['-O_ROUTES_FILE']
# if '-LOAD_ROUTES_FROM_FILE' in args:
#     if args['-LOAD_ROUTES_FROM_FILE'].upper() == 'TRUE':
#         LOAD_ROUTES_FROM_FILE = True
# if '-I_ROUTES_FILE' in args:
#     i_routes_file = args['-I_ROUTES_FILE']
# if '-MAX_SAMPLE_COUNT' in args:
#     max_sample_count = int(args['-MAX_SAMPLE_COUNT'])
# if '-SCORE_ROUTES' in args:
#     if args['-SCORE_ROUTES'].upper() == 'TRUE':
#         SCORE_ROUTES = True
# if '-I_GENE_DATA_FILE' in args:
#     i_gene_data_file = args['-I_GENE_DATA_FILE']
# if '-O_SCORES_FILE' in args:
#     o_scores_file = args['-O_SCORES_FILE']
# if '-O_FOLDER' in args:
#     o_folder = args['-O_FOLDER']
# if '-ZERO_NODE_PENALTY' in args:
#     zero_node_penalty = float(args['-ZERO_NODE_PENALTY'])
# if '-NODE_CONTRIBUTION_MAX' in args:
#     node_contribution_max = float(args['-NODE_CONTRIBUTION_MAX'])
# if '-GENE_DATA_DELIMITER' in args:
#     delimiter = args['-GENE_DATA_DELIMITER']
#     if delimiter.upper() == "COMMA":
#         gene_data_delimiter = ","
#     else:
#         gene_data_delimiter = "/t"

# if FIND_ROUTES:
#     o_routes_file_name = o_routes_file[:o_routes_file.rfind(".")]
#     o_routes_file_ext = o_routes_file[o_routes_file.rfind("."):]
#     if suffix != '' and suffix != 'NULL':
#         o_routes_file_name = o_routes_file_name + "_" + suffix    
#     o_routes_file = o_routes_file_name + o_routes_file_ext   

# if SCORE_ROUTES:
#     o_scores_file_name = o_scores_file[:o_scores_file.rfind(".")]
#     o_scores_file_ext = o_scores_file[o_scores_file.find("."):]
#     if suffix != '' and suffix != 'NULL':
#         o_scores_file_name = o_scores_file_name + "_" + suffix
#     o_scores_file =  o_scores_file_name + o_scores_file_ext   


i_pathways_folder = "C:/Users/whl19/Documents/Code/rPAC/pathways/"
o_routes_file = "../routes.txt"
i_routes_file=""

i_gene_data_file ="C:/Users/whl19/Documents/Code/BOSC2022/Dataset/Good/Bulkcell/combine_analysis/BulkcombinedLog2.txt" 
o_scores_file ="C:/Users/whl19/Documents/Code/BOSC2022/Dataset/Good/Bulkcell/combine_analysis/rpac_route_source_BulkcombinedLog2.txt" 

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

def main():
    
    if FIND_ROUTES:
        all_routes = route.find_routes_on_pathways(i_pathways_folder, route_end_candidate_types, True, o_routes_file)         
    elif LOAD_ROUTES_FROM_FILE:
        all_routes = route.load_routes_from_file(i_routes_file)
    
    if SCORE_ROUTES:
        score.set_node_contribution_max(node_contribution_max)
        score.set_zero_node_penalty(zero_node_penalty)
        
        score.score_routes_from_gene_data_file(i_pathways_folder, all_routes, i_gene_data_file, gene_data_delimiter,
                        max_sample_count, True, o_scores_file)
    '''
    if DRAW_ROUTES:
        visualize_routes(i_pathways_folder, all_routes, o_image_folder)
    '''

'''
def visualize_routes(pathways_folder, all_routes, o_image_folder):
    header_row = all_routes[0]
    for route in all_routes[1:]:
        route_id = route[0]
        if route_id != '8303':
            continue
        o_image_file = o_image_folder + "route" + route_id + ".png"
        rpax_pathway_graph.visualize_route(pathways_folder, route, header_row, o_image_file)
'''

if __name__ == '__main__':
    main()