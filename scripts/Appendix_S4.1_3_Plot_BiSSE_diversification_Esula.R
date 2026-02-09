# Plotting BiSSE Esula diversification rates
library(ggplot2)
source("Appendix_S4_multiplot.R")
data <- read.table("output_BiSSE_Esula_15000generaciones/BiSSE_Esula.log", header=TRUE)

start <- round(0.5*length(data$extinction.1))
end <- length(data$extinction.1)
BiSSE_types <- rep(c("1", "2"), each = length(data$extinction.1[start:
end]))
BiSSE_types2 <- rep(c("12", "21"), each = length(data$extinction.1[start:
end]))

# Extinction rate
dat_ext <- data.frame(dens = c(data$extinction.1[start:end], data$extinction.2[start:
end]), Type =
BiSSE_types)

# Speciation rate
dat_spec <- data.frame(dens = c(data$speciation.1[start:end], data$speciation.2[start:
end]), Type =
BiSSE_types)

# Dispersal rate
dat_disp<- data.frame(dens = c(data$rate_12[start:end], data$rate_21[start:
end]), Type =BiSSE_types2)

# Diversification rate
dat_div <- data.frame(dens = c(data$speciation.1[start:end]-data$extinction.1[start:
end], data$speciation.2[start:end]-data$extinction.2[start:end]), Type = BiSSE_types)

# Turnover
dat_rel <- data.frame(dens = c(data$extinction.1[start:end]/data$speciation.1[start:
end], data$extinction.2[start:end]/data$speciation.2[start:end]), Type = BiSSE_types)

pdf("Plot_BiSSE_diversification_Esula.pdf")
p1 <- ggplot(dat_spec, aes(x = dens, fill = Type)) + labs(title = "Speciation", x="
Rate", y="Posterior Density") + geom_density(alpha = 0.5) + guides(fill=
guide_legend(ncol=2,byrow=TRUE))
p2 <- ggplot(dat_ext, aes(x = dens, fill = Type)) + labs(title = "Extinction", x="Rate
", y="Posterior Density") + geom_density(alpha = 0.5) + guides(fill=guide_legend(
ncol=2,byrow=TRUE))
p3 <- ggplot(dat_div, aes(x = dens, fill = Type)) + labs(title = "Net-Diversification
", x="Rate", y="Posterior Density") + geom_density(alpha = 0.5) + guides(fill=
guide_legend(ncol=2,byrow=TRUE))
p4 <- ggplot(dat_rel, aes(x = dens, fill = Type)) + labs(title = "Relative Extinction
", x="Rate", y="Posterior Density") + geom_density(alpha = 0.5) + guides(fill=
guide_legend(ncol=2,byrow=TRUE))
p5 <- ggplot(dat_disp, aes(x = dens, fill = Type)) + labs(title = "Transition", x="
Rate", y="Posterior Density") + geom_density(alpha = 0.5) + guides(fill=
                                                                     guide_legend(ncol=2,byrow=TRUE))

multiplot(p1, p2, p3, p4, p5) 
dev.off()