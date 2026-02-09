### Plotting stochastic character mapping BiSSE ###
library(plotrix)
library(phytools)
library(ggtree)
library(ggplot2)
character_file = "output_BiSSE_Esula_15000generaciones/BiSSE_marginal_character_Esula.tree"

sim2 = read.simmap(file=character_file, format="phylip")
ladderize.simmap(sim2,right=FALSE)->sim2
plot(sim2)

# Define colours for the two states: annual and perennial
colors = vector()
for (i in 1:length( sim2$maps ) ) { 
    colors = c(colors, names(sim2$maps[[i]]) )
}
colors = sort(as.numeric(unique(colors)))
colors
cols = setNames( rainbow(length(colors), start=0.0, end=0.9), colors)
cols
# Change colours
cols[[1]] <- "#FFC20A"
cols[[2]] <- "#009E73"

# Plot stochastic character mapping results
# We do it three times, 1st without saving in PDF to get the legend, 2nd in PDF without labels to have the timeline adjusted to the X axis, 3rd in PDF with labels to copy and paste them in 2nd.
# 1st to get and copy the legend in the 2nd (it gives an error when adding it directly in the 2nd and 3rd)
plotSimmap(sim2, cols, fsize=0.2, direction="rightwards", lwd=1, split.vertical=TRUE, type="phylogram")
# Add legend
leg = names(cols)
leg
add.simmap.legend(leg, colors=cols, cex=0.3, x=0.8, y=0.8, fsize=0.8)
# Save the image using Plot > Export > Save as

# 2nd in PDF without labels to have the timeline adjusted to the X axis
pdf("Plot_BiSSE_stochastic_character_mapping_timeline_adjusted_no_labels_Esula.pdf", paper="special", height =16, width=11)
plotSimmap(sim2, cols, fsize=0, offset=NULL, direction="rightwards", lwd=1, split.vertical=TRUE, type="phylogram")
xaxp<-par()$xaxp[c(2,1,3)]
par(usr=par()$usr[c(2,1,3,4)])
par(xaxp=xaxp)
axis(1,cex.axis=0.8, pos = -1)
clip(x1=0,x2=42,y1=0,
     y2=Ntip(sim2))
dev.off()

# 3rd with labels to be copied and pasted into 2nd
pdf("Plot_BiSSE_stochastic_character_mapping_timeline_no_adjusted_labels_Esula.pdf", paper="special", height =16, width=11)
plotSimmap(sim2, cols, fsize=0.2, direction="rightwards", lwd=1, split.vertical=TRUE, type="phylogram")
xaxp<-par()$xaxp[c(2,1,3)]
par(usr=par()$usr[c(2,1,3,4)])
par(xaxp=xaxp)
axis(1,cex.axis=0.8, pos = -1)
clip(x1=0,x2=42,y1=0,
     y2=Ntip(sim2))
dev.off()



### Plotting posteriors of stochastic character mapping BiSSE ###
posterior_file = "output_BiSSE_Esula_15000generaciones/BiSSE_marginal_posterior_Esula.tree"

sim_p = read.simmap(file=posterior_file, format="phylip")
ladderize.simmap(sim_p,right=FALSE)->sim_p
plot(sim_p)

# Define colours for posterior probability 
colors = vector()
for (i in 1:length( sim_p$maps ) ) { 
  colors = c(colors, names(sim_p$maps[[i]]) )
}
colors = sort(as.numeric(unique(colors)))

# We can use different two colour choices to plot the posterior tree as a "heatmap". For posteriors, this works better
cols = setNames( heat.colors(length(colors), rev=TRUE), colors)

# Plot stochastic posterior mapping results
# We do it three times, 1st without saving in PDF to get the legend, 2nd in PDF without labels to have the timeline adjusted to the X axis, 3rd in PDF with labels to copy and paste them in 2nd.
# 1st to get and copy the legend in the 2nd (it gives error when adding it directly in the 2nd and 3rd)
plotSimmap(sim_p, colors=cols, fsize=0.2, direction="rightwards", lwd=1, split.vertical=TRUE, type="phylogram")
# Add legend
leg = names(cols)
leg
add.simmap.legend(leg, colors=cols, cex=0.3, x=0.2, y=0.2, fsize=0.3)
# Save the image using Plot > Export > Save as

# 2nd in PDF without labels to have the timeline adjusted to the X axis
pdf("Plot_BiSSE_stochastic_posterior_mapping_timeline_adjusted_no_labels_Esula.pdf", paper="special", height =16, width=11)
plotSimmap(sim_p, cols, fsize=0, offset=NULL, direction="rightwards", lwd=1, split.vertical=TRUE, type="phylogram")
xaxp<-par()$xaxp[c(2,1,3)]
par(usr=par()$usr[c(2,1,3,4)])
par(xaxp=xaxp)
axis(1,cex.axis=0.8, pos = -1)
clip(x1=0,x2=42,y1=0,
     y2=Ntip(sim_p))
dev.off()

# 3rd with labels to be copied and pasted into 2nd
pdf("Plot_BiSSE_stochastic_posterior_mapping_timeline_no_adjusted_labels_Esula.pdf", paper="special", height =16, width=11)
plotSimmap(sim_p, cols, fsize=0.2, direction="rightwards", lwd=1, split.vertical=TRUE, type="phylogram")
xaxp<-par()$xaxp[c(2,1,3)]
par(usr=par()$usr[c(2,1,3,4)])
par(xaxp=xaxp)
axis(1,cex.axis=0.8, pos = -1)
clip(x1=0,x2=42,y1=0,
     y2=Ntip(sim_p))
dev.off()