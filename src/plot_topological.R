setwd("/home/nat/Dropbox/cary_projects/DODA/cary/new_conv_2021/strict/recon_topo/br_part_1st_origin/topo") # your wd here
non_conv_siteLH <- read.csv("non_conv.sitelh.csv", header = F)
conv_siteLH <- read.csv("conv.sitelh.csv", header = F)

diffs <- conv_siteLH - non_conv_siteLH
scores <- exp(conv_siteLH - non_conv_siteLH)/(exp(conv_siteLH - non_conv_siteLH) + 1)
diffs$scores <- scores$V1
diffs$pos <- row.names(diffs)
colnames(diffs) <- c("deltaLogL", "score", "pos")
above <- diffs[which(diffs$score > 0.9),]
above_sorted <- above[order(above$score,decreasing = T),3]
above_str <- paste(above_sorted, collapse=", ")

plot(sort(diffs$deltaLogL,decreasing = T), col = ifelse(sort(diffs$deltaLogL,decreasing = T) < 0, "red", "blue"), type = "h",
     ylab=expression(paste(Delta," site logL, conv - non-conv")),
     xlab="", xaxt="n",
     main=expression(paste(Delta," Site Log-Likelihood, Convergent vs Non-Convergent Topology")),
     sub="PRANK, cleaned 0.1, JTT+G")  # change here 
abline(-2.2, 0, lty=3, col="red")
abline(2.2, 0, lty=3, col="blue")
text(length(diffs$deltaLogL)/2, 3.5, paste("Convergent (LR > 9, score > 0.9)\nSites:",above_str), col="blue")
text(length(diffs$deltaLogL)/2, -3.5, "Divergent (LR > 9, score > 0.9)", col="red")

