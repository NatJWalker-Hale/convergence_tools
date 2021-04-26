setwd("~/Dropbox/cary_projects/DODA/cary/new_conv_2021/strict/recon_topo/")

# diffsel

meandiffsel <- read.csv("diffsel/2cond/run1_cln0.1_1.meandiffsel", sep = "\t", header = F)[,1:21]
colnames(meandiffsel) <- c("pos",'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
score <- t(apply(meandiffsel, 1, function(x) {2 * abs(0.5 - x[2:21])}))
meandiffsel$max <- apply(score, 1, max)
meandiffsel$max_aa <- apply(score, 1, function(x) {colnames(score)[which.max(x)]})

# plotting

nsites <- max(meandiffsel$pos) # 273 in cleaned PRANK alignment
# 55 sites at a time, 55, 55, 55, 55, 53
par(mfrow=c(5,1), mar=c(2.0,0.5,0.5,0.5))
plot(meandiffsel$max[1:55], type="h", lwd=1, main="", xlab="", ylab="",
     xlim=c(1,55))
abline(0.9,0,lty=3)
points(which(meandiffsel$max[1:55] > 0.9), meandiffsel$max[meandiffsel$max > 0.9][1:55],
       pch=20)
plot(meandiffsel$max[56:110], type="h", lwd=1, main="", xlab="", ylab="",
     xaxt="n")
axis(1, at=seq(10,50,10), labels=c(65,75,85,95,105))
abline(0.9,0,lty=3)
points(which(meandiffsel$max[56:110] > 0.9), meandiffsel$max[meandiffsel$max[56:110] > 0.9],
       pch=20)
plot(meandiffsel$max[111:165], type="h", lwd=1, main="", xlab="", ylab="",
     xaxt="n")
axis(1, at=seq(10,50,10), labels=c(120,130,140,150,160))
abline(0.9,0,lty=3)
points(which(meandiffsel$max[111:165] > 0.9), meandiffsel$max[meandiffsel$max > 0.9][56:110],
       pch=20)
plot(meandiffsel$max[166:220], type="h", lwd=1, main="", xlab="", ylab="",
     xaxt="n")
axis(1, at=seq(10,50,10), labels=c(175,185,195,205,215))
abline(0.9,0,lty=3)
points(which(meandiffsel$max > 0.9), meandiffsel$max[meandiffsel$max > 0.9],
       pch=20)
plot(meandiffsel$max[221:273], type="h", lwd=1, main="", xlab="", ylab="",
     xaxt="n")
axis(1, at=seq(10,50,10), labels=c(230,240,250,260,270))
abline(0.9,0,lty=3)
points(which(meandiffsel$max > 0.9), meandiffsel$max[meandiffsel$max > 0.9],
       pch=20)

