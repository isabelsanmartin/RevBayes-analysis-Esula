library(diversitree)

tab <- read.table("output_HiSSE_Esula_15000generaciones/HiSSE_Esula.log", header=TRUE)
head(tab)

mcmc2diff <- with(tab, data.frame(s.diffA = speciation.1. - speciation.2., s.diffB = speciation.3. - speciation.4.,
                                  x.diffA = extinction.1. - extinction.2., x.diffB = extinction.3. - extinction.4., d.diff = rate_12 - rate_21))

head(mcmc2diff)

pdf("Plot_HiSSE_significant_differences_Esula.pdf", paper = "a4")
par(mfrow = c(2, 3), mar = c(3, 4, 2, 1)) # We calculate whether the estimated posterior probabilities for annual and perennial are significantly different

# Are these differences significant?
profiles.plot(mcmc2diff[1], col.line = "lightblue", xlab = "", ylab = "", main="s1A-s2A", xlim=c(-0.1,0.15))
abline(v=mean(as.numeric(as.character(mcmc2diff[,1]))), col="black", lwd=2)
abline(v=0, col="red", lwd=2)
profiles.plot(mcmc2diff[2], col.line = "lightblue", xlab = "", ylab = "", main="s1B-s2B", xlim=c(-0.5,0.5))
abline(v=mean(as.numeric(as.character(mcmc2diff[,2]))), col="black", lwd=2)
abline(v=0, col="red", lwd=2)
profiles.plot(mcmc2diff[3], col.line = "lightblue", xlab = "", ylab = "", main="x1A-x2A", xlim=c(-0.1,0.2))
abline(v=mean(as.numeric(as.character(mcmc2diff[,3]))), col="black", lwd=2)
abline(v=0, col="red", lwd=2)
profiles.plot(mcmc2diff[4], col.line = "lightblue", xlab = "", ylab = "", main="x1B-x2B", xlim=c(-0.1,0.2))
abline(v=mean(as.numeric(as.character(mcmc2diff[,4]))), col="black", lwd=2)
abline(v=0, col="red", lwd=2)
profiles.plot(mcmc2diff[5], col.line = "lightblue", xlab = "", ylab = "", main="d12-d21")
abline(v=mean(as.numeric(as.character(mcmc2diff[,5]))), col="black", lwd=2)
abline(v=0, col="red", lwd=2)
dev.off()