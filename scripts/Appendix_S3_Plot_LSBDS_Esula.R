library(RevGadgets)
library(treeio)
library(ggplot2)
setwd("/Users/imi/Downloads/9. RevBayes BDS Branch-Specific diversification rate estimation/BDS")
# Assign pointers to data files from RevBayes
treefile = "Esula_tree_no_outgroups.tre"
logfile = "BDS_Esula_rates.log"

# Read the tree
tree <- readTrees(paths = treefile)

# Read the log data with rates per branch
branch_data <- readTrace(logfile)

# Create annotated tree with the branch data
annotated_tree <- processBranchData(tree, branch_data, burnin = 0.25, parnames = c("avg_lambda", "avg_mu", "num_shifts"), summary = "median", net_div = TRUE)

p <- plotTree(tree = annotated_tree,
              node_age_bars = FALSE,
              node_pp = FALSE,
              tip_labels = TRUE,
              color_branch_by = "avg_lambda",
              line_width = 0.8,
              tip_labels_size = 1,
              branch_color = c("blue","green")) +
  ggplot2::theme(legend.position=c(.1, .9));p
ggsave("LSBDS_speciation_avg_lambda_Esula.pdf", width=40, height=49, units="cm")

p <- plotTree(tree = annotated_tree,
              node_age_bars = FALSE,
              node_pp = FALSE,
              tip_labels = TRUE,
              color_branch_by = "avg_mu",
              line_width = 0.8,
              tip_labels_size = 1,
              branch_color = c("blue","green")) +
     ggplot2::theme(legend.position=c(.1, .9));p
ggsave("LSBDS_extinction_avg_mu_Esula.pdf", width=40, height=49, units="cm") 