require(yaml)
require(viridis)
require(RColorBrewer)
require(grid)

setwd("~/Dropbox/cary_projects/DODA/cary/new_conv_2021/strict/recon_topo/") # set wd here

# diffsel

meandiffsel <- read.csv("br_part_special/diffsel/2cond/run1_cln_1.meandiffsel", sep = "\t", header = F)[,1:21]
colnames(meandiffsel) <- c("pos",'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
score <- t(apply(meandiffsel, 1, function(x) {2 * abs(0.5 - x[2:21])}))
meandiffsel$max <- apply(score, 1, max)
meandiffsel$max_aa <- apply(score, 1, function(x) {colnames(score)[which.max(x)]})

# msd

msd <- read.csv("br_part_original/prank/msd/output.txt", sep="\t", skip=3, header=F)
msd <- msd[order(msd$V1),]
msd$score <- 1-as.numeric(msd$V2)

# PCOC

pcoc <- read.csv("br_part_special/PCOC/prank_PCOC_C60_aa_brlen_unroot/RUN_20210430_084954/DODAa_combined_no_og_strict_for_synth.cds.fa.nostop.name.noF.best.fas.results.tsv",
                 sep="\t", header=T) # no transformations necessary

# tdg09

out <- yaml.load_file("br_part_special/prank/tdg09/prank_JTT+G_out.txt")
tdg09 <- as.data.frame(matrix(unlist(out$FullResults), ncol=9, byrow=T))
colnames(tdg09) <- c("pos","ssFparam","ssFlogL","lssFparam","lssFlogL","deltaLogL","df","LRT","FDR")
tdg09$score <- 1-as.numeric(tdg09$FDR)
tdg09$score[is.na(tdg09$score)] <- 0.0

# topo

topo <- read.csv("br_part_special/prank/topo/prank_JTT+G_topo_out.txt", header=T)


# all

results <- data.frame(pos=pcoc$Sites,
                      diffsel=meandiffsel$max,
                      msd=msd$score,
                      pcoc=pcoc$PCOC,
                      tdg09=tdg09$score,
                      topo=topo$topological
                      )

results$pass <- apply(results[,2:6],1,function(x) {sum(x > 0.9)})

#write.csv(x = results, "br_part_original/prank/all_results.csv", row.names = F)
#write.csv(x = results, "br_part_1st_origin/prank/all_results.csv", row.names = F)
#write.csv(x = results, "br_part_special/prank/all_results.csv", row.names = F)
#write.csv(x = results,"br_part_original_inverse/prank/all_results.csv", row.names = F)

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

require(VennDiagram)

# prank
df <- read.csv("br_part_original/prank/all_results.csv", header=T)
original_pass3 <- df[which(df$pass >= 3),]$pos
df <- read.csv("br_part_1st_origin/prank/all_results.csv", header=T)
first_origin_pass3 <- df[which(df$pass >= 3),]$pos
df <- read.csv("br_part_special/prank/all_results.csv", header=T)
special_pass3 <- df[which(df$pass >= 3),]$pos
df <- read.csv("br_part_original_inverse/prank/all_results.csv", header=T)
inverse_pass3 <- df[which(df$pass >= 3),]$pos

sites <- list(Original=original_pass3, First=first_origin_pass3, Special=special_pass3, Inverse=inverse_pass3)

venn.diagram(sites, filename = "test.png")

display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(sites,
             fill = brewer.pal(4, "Set2"),
             lwd=2,
             lty="blank",
             sub.fontface="arial",
             sub.fontfamily="sanserif")

inAll <- intersect(intersect(intersect(original_pass3,first_origin_pass3),special_pass3),inverse_pass3)
uniq_original <- c(78)
uniq_special <- c(NA)
uniq_inverse <- c(183)
uniq_first <- c(79,81,148)
