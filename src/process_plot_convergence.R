require(yaml)
require(viridis)
require(RColorBrewer)
require(grid)

setwd("~/Dropbox/cary_projects/DODA/cary/new_conv_2021/strict/recon_topo/") # set wd here

# diffsel

meandiffsel <- read.csv("br_part_1st_origin/diffsel/2cond/run1_cln_1.meandiffsel", sep = "\t", header = F)[,1:21]
colnames(meandiffsel) <- c("pos",'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
score <- t(apply(meandiffsel, 1, function(x) {2 * abs(0.5 - x[2:21])}))
meandiffsel$max <- apply(score, 1, max)
meandiffsel$max_aa <- apply(score, 1, function(x) {colnames(score)[which.max(x)]})

# msd

msd <- read.csv("br_part_original/msd/output.txt", sep="\t", skip=3, header=F)
msd <- msd[order(msd$V1),]
msd$score <- 1-as.numeric(msd$V2)

# PCOC

pcoc <- read.csv("br_part_1st_origin/PCOC/prank_PCOC_C60_aa_brlen_unroot_exclude/RUN_20210429_094808/DODAa_combined_no_og_strict_for_synth.cds.fa.nostop.name.noF.best.fas.results.tsv",
                 sep="\t", header=T) # no transformations necessary

# tdg09

out <- yaml.load_file("br_part_1st_origin/tdg09/tdg09_prank_JTT+G_exclude_out.txt")
tdg09 <- as.data.frame(matrix(unlist(out$FullResults), ncol=9, byrow=T))
colnames(tdg09) <- c("pos","ssFparam","ssFlogL","lssFparam","lssFlogL","deltaLogL","df","LRT","FDR")
tdg09$score <- 1-as.numeric(tdg09$FDR)
tdg09$score[is.na(tdg09$score)] <- 0.0

# topo

topo <- read.csv("br_part_1st_origin/topo/prank_topo_out.txt", header=T)


# all

results <- data.frame(pos=pcoc$Sites,
                      diffsel=meandiffsel$max,
                      msd=msd$score,
                      pcoc=pcoc$PCOC,
                      tdg09=tdg09$score,
                      topo=topo$topological
                      )

results$pass <- apply(results[,2:6],1,function(x) {sum(x > 0.9)})

cols <- brewer.pal(5, "Set2")
# cols <- viridis(5)

# plotting

sitedf <- data.frame(method=c("diffsel","msd","pcoc","tdg09","topo"),
                      score=as.numeric(results[22,][2:6]))
par(bg="lightgrey")
barplot(sitedf$score, axes=F, ylim=c(0,1), col = cols)
mtext(text = 2, side = 1, line = 1)

# pdf("test.pdf", paper="a4")
par(mfrow=c(10,28), mar=c(2.0,0.5,0.5,0.5))
for (i in 1:max(results$pos)) {
  sitedf <- data.frame(method=c("diffsel","msd","pcoc","tdg09","topo"),
                       score=as.numeric(results[i,][2:6]))
  # if (sum(sitedf$score > 0.9) == 4 ) {
  #   par(bg="gold")
  # }
  if (i %in% c(1,29,57,85,113,141,169,197,225,253)) {
    barplot(sitedf$score, ylim=c(0,1), col = cols)
  } else {
  barplot(sitedf$score, yaxt="n", ylim=c(0,1), col = cols)
  }
  # if (i %% 5 == 0) {
    mtext(text=i, side=1, line = 1, cex = 0.7)
  # }
  # if (i %in% c(28,56,84,112,140,168,196,224,252,274)) {
  #   (Y <- grconvertY(0.8, "user", "nic"))
  #   pushViewport(viewport())
  #   grid.lines(x = c(0,4), y = Y, gp = gpar(lty=3, col = "black"))
  #   popViewport()
  # }
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
#legend("bottom",legend=c("diffsel","msd","pcoc","tdg09","topo"), col = cols)
legend('bottom',legend=c("diffsel","msd","pcoc","tdg09","topo"), col=cols, lwd = 5, xpd = TRUE, horiz = TRUE, cex = 1, seg.len=1, bty = 'n')


(Y <- grconvertY(0.9, "user", "ndc"))
pushViewport(viewport())
grid.lines(x = c(0,4), y = Y, gp = gpar(lty=3, col = "black"))
popViewport()

# nsites <- max(meandiffsel$pos) # 274 in cleaned PRANK alignment
# # 55 sites at a time, 55, 55, 55, 55, 55
# par(mfrow=c(5,1), mar=c(2.0,0.5,0.5,0.5))
# plot(meandiffsel$max[1:55], type="h", lwd=1, main="", xlab="", ylab="",
#      xlim=c(1,55))
# abline(0.9,0,lty=3)
# # points(which(meandiffsel$max[1:55] > 0.9), meandiffsel$max[meandiffsel$max > 0.9][1:55],
# #        pch=20)
# plot(meandiffsel$max[56:110], type="h", lwd=1, main="", xlab="", ylab="",
#      xaxt="n")
# axis(1, at=seq(10,50,10), labels=c(65,75,85,95,105))
# abline(0.9,0,lty=3)
# # points(which(meandiffsel$max[56:110] > 0.9), meandiffsel$max[meandiffsel$max[56:110] > 0.9],
# #        pch=20)
# plot(meandiffsel$max[111:165], type="h", lwd=1, main="", xlab="", ylab="",
#      xaxt="n")
# axis(1, at=seq(10,50,10), labels=c(120,130,140,150,160))
# abline(0.9,0,lty=3)
# # points(which(meandiffsel$max[111:165] > 0.9), meandiffsel$max[meandiffsel$max > 0.9][56:110],
# #        pch=20)
# plot(meandiffsel$max[166:220], type="h", lwd=1, main="", xlab="", ylab="",
#      xaxt="n")
# axis(1, at=seq(10,50,10), labels=c(175,185,195,205,215))
# abline(0.9,0,lty=3)
# # points(which(meandiffsel$max > 0.9), meandiffsel$max[meandiffsel$max > 0.9],
# #        pch=20)
# plot(meandiffsel$max[221:273], type="h", lwd=1, main="", xlab="", ylab="",
#      xaxt="n")
# axis(1, at=seq(10,50,10), labels=c(230,240,250,260,270))
# abline(0.9,0,lty=3)
# # points(which(meandiffsel$max > 0.9), meandiffsel$max[meandiffsel$max > 0.9],
# #        pch=20)

