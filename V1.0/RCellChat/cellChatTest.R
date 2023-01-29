# Load the required libraries

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

###################################################################################################################################################
# Part I: Data input & processing and initialization of CellChat object
# CellChat requires two user inputs: one is the gene expression data of cells, and the other is either 
# user assigned cell labels (i.e., label-based mode) or a low-dimensional representation of the single-cell data 
# (i.e., label-free mode).
# For the latter, CellChat automatically groups cells by building a shared neighbor graph based on the cell-cell 
# distance in the low-dimensional space or the pseudotemporal trajectory space.
###########################################################################################################################################################################
# Load data
# For the gene expression data matrix, genes should be in rows with rownames and cells in columns with colnames. 
# Normalized data (e.g., library-size normalization and then log-transformed with a pseudocount of 1) is required as input for CellChat analysis. 
# If user provides count data, we provide a normalizeData function to account for library size and then do log-transformed. 
# For the cell group information, a dataframe with rownames is required as input for CellChat.

# Here we load a scRNA-seq data matrix and its associated cell meta data
load(url("https://ndownloader.figshare.com/files/25950872")) # This is a combined data from two biological conditions: normal and diseases
# load("/Users/suoqinjin/Documents/CellChat/tutorial/data_humanSkin_CellChat.rda")
data.input = data_humanSkin$data # normalized data matrix
meta = data_humanSkin$meta # a dataframe with rownames containing cell mata data
cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data

# Prepare input data for CelChat analysis
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
unique(meta$labels) # check the cell labels
# ###########################################################################################################################################################################
# Create a CellChat object
# USERS can create a new CellChat object from a data matrix, Seurat or SingleCellExperiment object. 
# If input is a Seurat or SingleCellExperiment object, the meta data in the object will be used by default and USER must provide group.
# by to define the cell groups. e.g, group.by = "ident" for the default cell identities in Seurat object.
# 
# NB: If USERS load previously calculated CellChat object (version < 0.5.0), please update the object via updateCellChat
##########################################################################################################################################################################
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")


##########################################################################################################################################################################
# Add cell information into meta slot of the object (Optional)
# If cell mata information is not added when creating CellChat object, 
# USERS can also add it later using addMeta, and set the default cell identities using setIdent.
##########################################################################################################################################################################
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

##########################################################################################################################################################################

# Set the ligand-receptor interaction database
# Our database CellChatDB is a manually curated database of literature-supported ligand-receptor interactions in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60% of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions, including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions and 16.5% of cell-cell contact interactions.
# 
# Users can update CellChatDB by adding their own curated ligand-receptor pairs.Please check our tutorial on how to do it.
##########################################################################################################################################################################
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

##########################################################################################################################################################################

# Preprocessing the expression data for cell-cell communication analysis
# To infer the cell state-specific communications, we identify over-expressed ligands or receptors in one cell group and then identify over-expressed ligand-receptor interactions if either ligand or receptor is over-expressed.
# 
# We also provide a function to project gene expression data onto protein-protein interaction (PPI) network. 
# Specifically, a diffusion process is used to smooth genes' expression values based on their neighbors' defined in a high-confidence experimentally validated protein-protein network.
# This function is useful when analyzing single-cell data with shallow sequencing depth because the projection reduces the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors. 
# One might be concerned about the possible artifact introduced by this diffusion process, however, it will only introduce very weak communications. 
# USERS can also skip this step and set raw.use = TRUE in the function computeCommunProb().
##########################################################################################################################################################################

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
 future::plan("multiprocess", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

##########################################################################################################################################################################

# Part II: Inference of cell-cell communication network
# CellChat infers the biologically significant cell-cell communication by assigning each interaction with a probability value and peforming a permutation test. 
# CellChat models the probability of cell-cell communication by integrating gene expression with prior known knowledge of the interactions between signaling ligands, receptors and their cofactors using the law of mass action.
# 
# The number of inferred ligand-receptor pairs clearly depends on the method for calculating the average gene expression per cell group. 
# By default, CellChat uses a statistically robust mean method called 'trimean', which produces fewer interactions than other methods. 
# However, we find that CellChat performs well at predicting stronger interactions, which is very helpful for narrowing down on interactions for further experimental validations.
# In computeCommunProb, we provide an option for using other methods, such as 5% and 10% truncated mean, to calculating the average gene expression. 
# Of note, 'trimean' approximates 25% truncated mean, implying that the average gene expression is zero if the percent of expressed cells in one group is less than 25%. 
# To use 10% truncated mean, USER can set type = "truncatedMean" and trim = 0.1. 
# If very well-known signaling pathways in the studied biological process are not predicted, USER can try truncatedMean with different trim values. 
# The function computeAveExpr can help to check the average expression of signaling genes of interest, e.g, computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1).
# 
# When analyzing unsorted single-cell transcriptomes, under the assumption that abundant cell populations tend to send collectively stronger signals than the rare cell populations, 
# CellChat can also consider the effect of cell proportion in each cell group in the probability calculation.
# USER can set population.size = TRUE.
# 
# Compute the communication probability and infer cellular communication network
##########################################################################################################################################################################
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
##########################################################################################################################################################################

# Extract the inferred cellular communication network as a data frame
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# 
# df.net <- subsetCommunication(cellchat) returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# 
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# 
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
# 
# Infer the cell-cell communication at a signaling pathway level
# CellChat computes the communication probability on signaling pathway level by summarizing the communication probabilities of all ligands-receptors interactions associated with each signaling pathway.
# 
# NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot 'net' and 'netP', respectively.
##########################################################################################################################################################################
cellchat <- computeCommunProbPathway(cellchat)
##########################################################################################################################################################################

# Calculate the aggregated cell-cell communication network
# We can calculate the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability. USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.
##########################################################################################################################################################################
cellchat <- aggregateNet(cellchat)
##########################################################################################################################################################################

# We can also visualize the aggregated cell-cell communication network. 
# For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
##########################################################################################################################################################################
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
##########################################################################################################################################################################

# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.
##########################################################################################################################################################################
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

##########################################################################################################################################################################

# Part III: Visualization of cell-cell communication network
# Upon infering the cell-cell communication network, CellChat provides various functionality for further data exploration, analysis, and visualization.
# 
# It provides several ways for visualizing cell-cell communication network, including hierarchical plot, circle plot, Chord diagram, and bubble plot.
# 
# It provides an easy-to-use tool for extracting and visualizing high-order information of the inferred networks. For example, it allows ready prediction of major signaling inputs and outputs for cell populations and how these populations and signals coordinate together for functions.
# 
# It can quantitatively characterize and compare the inferred cell-cell communication networks using an integrated approach by combining social network analysis, pattern recognition, and manifold learning approaches.
# 
# Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
# Hierarchy plot: USER should define vertex.receiver, which is a numeric vector giving the index of the cell groups as targets in the left part of hierarchy plot. This hierarchical plot consist of two components: the left portion shows autocrine and paracrine signaling to certain cell groups of interest (i.e, the defined vertex.receiver), and the right portion shows autocrine and paracrine signaling to the remaining cell groups in the dataset. Thus, hierarchy plot provides an informative and intuitive way to visualize autocrine and paracrine signaling communications between cell groups of interest. For example, when studying the cell-cell communication between fibroblasts and immune cells, USER can define vertex.receiver as all fibroblast cell groups.
# 
# Chord diagram: CellChat provides two functions netVisual_chord_cell and netVisual_chord_gene for visualizing cell-cell communication with different purposes and different levels. netVisual_chord_cell is used for visualizing the cell-cell communication between different cell groups (where each sector in the chord diagram is a cell group), and netVisual_chord_gene is used for visualizing the cell-cell communication mediated by mutiple ligand-receptors or signaling pathways (where each sector in the chord diagram is a ligand, receptor or signaling pathway.)
# 
# Explnations of edge color/weight, node color/size/shape: In all visualization plots, edge colors are consistent with the sources as sender, and edge weights are proportional to the interaction strength. Thicker edge line indicates a stronger signal. In the Hierarchy plot and Circle plot, circle sizes are proportional to the number of cells in each cell group. In the hierarchy plot, solid and open circles represent source and target, respectively. In the Chord diagram, the inner thinner bar colors represent the targets that receive signal from the corresponding outer bar. The inner bar size is proportional to the signal strength received by the targets. Such inner bar is helpful for interpreting the complex chord diagram. Note that there exist some inner bars without any chord for some cell groups, please just igore it because this is an issue that has not been addressed by circlize package.
# 
# Visualization of cell-cell communication at different levels: One can visualize the inferred communication network of signaling pathways using netVisual_aggregate, and visualize the inferred communication networks of individual L-R pairs associated with that signaling pathway using netVisual_individual.
# 
# Here we take input of one signaling pathway as an example. All the signaling pathways showing significant communications can be accessed by cellchat@netP$pathways.
##########################################################################################################################################################################

pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

##########################################################################################################################################################################
# For the chord diagram, CellChat has an independent function netVisual_chord_cell to flexibly visualize the signaling network by adjusting different parameters in the circlize package. For example, we can define a named char vector group to create multiple-group chord diagram, e.g., grouping cell clusters into different cell types.
##########################################################################################################################################################################
# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

##########################################################################################################################################################################
# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
##########################################################################################################################################################################



netAnalysis_contribution(cellchat, signaling = pathways.show)









