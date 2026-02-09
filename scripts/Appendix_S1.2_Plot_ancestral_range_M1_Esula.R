setwd("C:/Users/Bird_/Downloads/Esula_Manuscrito/New_Analysis/Stratified_M1")
install.packages("devtools", dependencies=TRUE)
install.packages("glue")
library(devtools)
install_github("GuangchuangYu/ggtree")
install_github("revbayes/RevGadgets")
library(RevGadgets)
library(ggtree)
source("Appendix_S1_Plot_ancestral_range_useful.R")
getwd()

# File names
fp = "./" # Edit to provide an absolute filepath
plot_fn = paste(fp, "output_esula_stratified_M1/stratified_M1.range.pdf",sep="")
tree_fn = paste(fp, "output_esula_stratified_M1/anc_states.ase.tre", sep="")
label_fn = paste(fp, "output_esula_stratified_M1/output_10000gen_epoch-nodist_ML_stratified_M1.state_labels.txt", sep="")

sessionInfo()

# Get state labels and state colours
labs <- c("1" = "AFR"
          ,"2" = "WP"
          ,"3" = "EP"
          ,"4" = "MAC"
          ,"5" = "NAM"
          ,"6" = "AFR-WP"
          ,"7" = "AFR-EP"
          ,"8" = "WP-EP"
          ,"9" = "AFR-MAC"
          ,"10" = "WP-MAC"
          ,"11" = "EP-MAC"
          ,"12" = "AFR-NAM"
          ,"13" = "WP-NAM"
          ,"14" = "EP-NAM"
          ,"15" = "MAC-NAM"
          ,"16" = "AFR-WP-EP"
          ,"17" = "AFR-WP-MAC"
          ,"18" = "AFR-EP-MAC"
          ,"19" = "WP-EP-MAC"
          ,"20" = "AFR-WP-NAM"
          ,"21" = "AFR-EP-NAM"
          ,"22" = "WP-EP-NAM"
          ,"23" = "AFR-MAC-NAM"
          ,"24" = "WP-MAC-NAM"
          ,"25" = "EP-MAC-NAM"
          ,"26" = "AFR-WP-EP-MAC"
          ,"27" = "AFR-WP-EP-NAM"
          ,"28" = "AFR-WP-MAC-NAM"
          ,"29" = "AFR-EP-MAC-NAM"
          ,"30" = "WP-EP-MAC-NAM"
          ,"31" = "AFR-WP-EP-MAC-NAM")

ancstates <- processAncStates(tree_fn, state_labels = labs)

# Plot the ancestral states
pp=plotAncStatesPie(t = ancstates, 
                    # Include cladogenetic events
                    cladogenetic = T,
                    # Add text labels to the tip pie symbols
                    tip_labels_states = FALSE,
                    # tama?o del nombre de las sp.
                    tip_labels_size = 4,
                    # Offset those text labels slightly
                    tip_labels_states_offset = .05,
                    # Offset the tip labels to make room for tip pies
                    tip_labels_offset = 1, 
                    # Move tip pies right slightly 
                    tip_pie_nudge_x = .07,
                    # Change the size of node and tip pies  
                    tip_pie_size = 0.25,
                    pie_colors = c("#d95f02","#F94C10","#B9B4C7","#016A70","#00FF7F","#3457D5","#D4ECDD","#C5A880","#FF7272","#a6d85f","#750050","#750050","#e7298a","#00FFFF","#783F04","#F4CCCC","#FFD523","#512B81","#CE2029","#394867","#99627A"),          
                    node_pie_size = 0.25,
                    state_transparency = 1)

# Get plot dimensions
x_phy = max(pp$data$x) # Get height of tree
x_label = 3.5 # Choose space for tip labels
x_start = 45 # Choose starting age (greater than x_phy)
x0 = -(x_start - x_phy) # Determine starting pos for xlim
x1 = x_phy + x_label # Determine ending pos for xlim

# Add axis
pp = pp + theme_tree2()
pp = pp + labs(x="Age (Ma)")

# Change x coordinates
pp = pp + coord_cartesian(xlim=c(x0,x1), expand=TRUE)

# Plot axis ticks, the values of 35, 25 and 3.5 are the millions of years, you can put whatever you want, just like the -A, +T, +E
island_axis = sec_axis(~ ., breaks=x_phy-c(35, 25, 3.5), labels=c("-A","+T","+E") )
x_breaks = seq(0,x_start,5) + x0
x_labels = rev(seq(0,x_start,5))
pp = pp + scale_x_continuous(breaks=x_breaks, labels=x_labels, sec.axis=island_axis)
pp

# Save, pie charts are resized with height and widht
ggsave(file=plot_fn, plot=pp, device="pdf", height=49, width=40, useDingbats=F)