# Read in data
setwd("~/Documents/HT_manuscript/divergencePlots/")
test <- read.table("Bos.taurus_bothL1andBovB.fasta.divsum", header=T)

# Draw boxplot
#mat <- matrix(data = NA, nrow = 10000,ncol = 10)
#mat[,ncol(mat)] <- sample(x = test$Div, size = 1000, replace = TRUE, prob = test$L1)
#mat[,ncol(mat)-5] <- sample(x = test$Div, size = 1000, replace = TRUE, prob = test$Tx1)
#mat[,ncol(mat)-2] <- sample(x = test$Div, size = 1000, replace = TRUE, prob = test$BovB)

#boxplot(mat, horizontal = TRUE, ylim = c(0,60), yaxs = "i", col = c("darkblue", "lightblue", "orange"),axes = FALSE, notch = TRUE)
#par(new = TRUE)

# Plot divergence bars
# Change the position of each bar by adding test$Div+0.25 or -0.25 to plot(_)
pdf("Bos.taurus_bothL1andBovB.pdf", width=6,height=5)
plot(test$L1, type="h", col="darkblue",lwd=4, 
     xlab="Kimura substitution level (CpG adjusted)",
     ylab="Coverage",
     main="Bos taurus", ylim=c(0,max(test$BovB)), xaxs = "i",xlim = c(-0.5,60))
points(test$Div+0.4, test$Tx12, type="h", col="lightblue", lwd=4)
points(test$Div+0.7, test$BovB, type="h", col="orange", lwd=4)
legend(47,max(test$BovB), c("L1", "Tx1", "BovB"),fill=c("darkblue", "lightblue", "orange"), box.lwd = 0,box.col = "white", bg="white", cex=1)
dev.off()

