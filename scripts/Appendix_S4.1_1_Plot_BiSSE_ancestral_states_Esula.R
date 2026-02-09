library(RevGadgets)
library(ape)

tree_file <- paste0("output_BiSSE_Esula_15000generaciones/anc_states_BiSSE_Esula.tree")
tree_file

# Process file
example <- processAncStates(tree_file,
                            state_labels = c("0" = "Annual",
                                             "1" = "Perennial"))

# Have states vary by color and indicate state pp with size (default)
plot <- plotAncStatesMAP(example)
plot