require(yaml)
require(viridis)
require(RColorBrewer)
require(grid)

setwd("/home/nat/Dropbox/cary_projects/DODA/new_conv_2021/strict/recon_topo/prank/brlen_cln_aln/br_part_1st_origin/exclude_kewa/") # set wd here

# diffsel

meandiffsel <- read.csv("diffsel/run1_1.meandiffsel", sep = "\t", header = F)[,1:21]
colnames(meandiffsel) <- c("pos",'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
score <- t(apply(meandiffsel, 1, function(x) {2 * abs(0.5 - x[2:21])}))
meandiffsel$max <- apply(score, 1, max)
meandiffsel$max_aa <- apply(score, 1, function(x) {colnames(score)[which.max(x)]})

# msd

msd <- read.csv("msd/output.txt", sep="\t", skip=3, header=F)
msd <- msd[order(msd$V1),]
msd$score <- 1-as.numeric(msd$V2)

# PCOC

pcoc <- read.csv(Sys.glob("PCOC_new1/RUN_*/*_results.tsv"),
                 sep="\t", header=T)
# switch NAs for too gappy sites to 0

pcoc$PCOC_V1[is.na(pcoc$PCOC_V1)] <- 0.0

# tdg09

out <- yaml.load_file("tdg09/tdg09_out.txt")
tdg09 <- as.data.frame(matrix(unlist(out$FullResults), ncol=9, byrow=T))
colnames(tdg09) <- c("pos","ssFparam","ssFlogL","lssFparam","lssFlogL","deltaLogL","df","LRT","FDR")
tdg09$score <- 1-as.numeric(tdg09$FDR)
tdg09$score[is.na(tdg09$score)] <- 0.0

# PELICAN

pelican <- read.table("PELICAN/all_sites.tsv", header=T)
pelican$score <- 1 - as.numeric(pelican$aagtr_pval)
pelican$score

# Godon

godon <- read.table("~/Dropbox/cary_projects/DODA/20220215_selection/strict/recon_topo/prank/all_foreground/BEB_results.tsv",
                    header=T)

# topo

topo <- read.csv("topo/topo_output.txt", header=T)


# all

results <- data.frame(pos=pcoc$Sites,
                      diffsel=meandiffsel$max,
                      msd=msd$score,
                      pcoc=pcoc$PCOC,
                      tdg09=tdg09$score,
                      topo=topo$topological
                      )

results <- data.frame(pos=pelican$site,
                      godon=godon$p[match(pelican$site, godon$pos)],
                      pcoc=pcoc$PCOC[match(pelican$site, pcoc$Sites)],
                      pelican=pelican$score,
                      diffsel=meandiffsel$max
                      )

results$godon[is.na(results$godon)] <- 0.0
results$pcoc[is.na(results$pcoc)] <- 0.0
 
write.csv(results, "all_output.csv", row.names = FALSE)

results$pass <- apply(results[,2:5],1,function(x) {sum(x > 0.95)})

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
    barplot(sitedf$score, ylim=c(0,1), col = cols, border = NA)
  } else {
  barplot(sitedf$score, yaxt="n", ylim=c(0,1), col = cols, border = NA)
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

# multiple separate bars on one line
# par(mfrow=c(5,1))

conv_sub_sites <- read.table("~/Dropbox/cary_projects/DODA/figures/thesis_versions/element_files/convergent_sub_sites_cln.txt", header=T)

conv_sub_sites_vec <- conv_sub_sites$pos_cln

col_vec <- rep("#cccccc", length(results$pos))
col_vec[conv_sub_sites_vec] <- conv_sub_sites$col

# diffsel

plot(results$diffsel, type="h", col=col_vec,
     frame.plot = F, xlim=c(1,274), ylim=c(0,1.0),
     ylab="", xlab="", xaxt="n")

pos0.95 <- results$pos[which(results$diffsel >= 0.95)]
points_col_vec <- rep("#cccccc", length(pos0.95))

for (i in pos0.95) {
  if (i %in% conv_sub_sites$pos_cln) {
    points_col_vec[which(pos0.95 == i)] <- conv_sub_sites$col[which(conv_sub_sites$pos_cln == i)]
  }
}

points(which(results$diffsel >= 0.95),
       results$diffsel[results$diffsel >= 0.95],
       col="#cccccc",
       pch=20)

points(pos0.95,
       results[which(results$pos %in% pos0.95),]$diffsel,
       col=points_col_vec,
       pch=20)

abline(h=0.95,
       lty=3)

# pelican

plot(results$pelican, type="h", col=col_vec,
     frame.plot = F, xlim=c(1,274), ylim=c(0,1.0),
     ylab="", xlab="", xaxt="n")

pos0.95 <- results$pos[which(results$pelican >= 0.95)]
points_col_vec <- rep("#cccccc", length(pos0.95))

for (i in pos0.95) {
  if (i %in% conv_sub_sites$pos_cln) {
    points_col_vec[which(pos0.95 == i)] <- conv_sub_sites$col[which(conv_sub_sites$pos_cln == i)]
  }
}

points(which(results$pelican >= 0.95),
       results$pelican[results$pelican >= 0.95],
       col="#cccccc",
       pch=20)

points(pos0.95,
       results[which(results$pos %in% pos0.95),]$pelican,
       col=points_col_vec,
       pch=20)

abline(h=0.95,
       lty=3)

# godon

plot(results$godon, type="h", col=col_vec,
     frame.plot = F, xlim=c(1,274), ylim=c(0,1.0),
     ylab="", xlab="", xaxt="n")

pos0.95 <- results$pos[which(results$godon >= 0.95)]
points_col_vec <- rep("#cccccc", length(pos0.95))

for (i in pos0.95) {
  if (i %in% conv_sub_sites$pos_cln) {
    points_col_vec[which(pos0.95 == i)] <- conv_sub_sites$col[which(conv_sub_sites$pos_cln == i)]
  }
}

points(which(results$godon >= 0.95),
       results$godon[results$godon >= 0.95],
       col="#cccccc",
       pch=20)

points(pos0.95,
       results[which(results$pos %in% pos0.95),]$godon,
       col=points_col_vec,
       pch=20)

abline(h=0.95,
       lty=3)

# msd

plot(results$msd, type="h", col=col_vec,
     frame.plot = F, xlim=c(1,274), ylim=c(0,1.0),
     ylab="", xlab="", xaxt="n")

pos0.95 <- results$pos[which(results$msd >= 0.95)]
points_col_vec <- rep("#cccccc", length(pos0.95))

for (i in pos0.95) {
  if (i %in% conv_sub_sites$pos_cln) {
    points_col_vec[which(pos0.95 == i)] <- conv_sub_sites$col[which(conv_sub_sites$pos_cln == i)]
  }
}

points(which(results$msd >= 0.95),
       results$msd[results$msd >= 0.95],
       col="#cccccc",
       pch=20)

points(pos0.95,
       results[which(results$pos %in% pos0.95),]$msd,
       col=points_col_vec,
       pch=20)

abline(h=0.95,
       lty=3)

# pcoc

plot(results$pcoc, type="h", col=col_vec,
     frame.plot = F, xlim=c(1,274), ylim=c(0,1.0),
     ylab="", xlab="", xaxt="n", yaxt="n")
axis(2, c(0, 0.2, 0.4, 0.6, 0.8, 1.0))

pos0.95 <- results$pos[which(results$pcoc >= 0.95)]
points_col_vec <- rep("#cccccc", length(pos0.95))

for (i in pos0.95) {
  if (i %in% conv_sub_sites$pos_cln) {
    points_col_vec[which(pos0.95 == i)] <- conv_sub_sites$col[which(conv_sub_sites$pos_cln == i)]
  }
}

points(which(results$pcoc >= 0.95),
       results$pcoc[results$pcoc >= 0.95],
       col="#cccccc",
       pch=20)

points(pos0.95,
       results[which(results$pos %in% pos0.95),]$pcoc,
       col=points_col_vec,
       pch=20)

abline(h=0.95,
       lty=3)

# tdg09

plot(results$tdg09, type="h", col=col_vec,
     frame.plot = F, xlim=c(1,274), ylim=c(0,1.0),
     ylab="", xlab="", xaxt="n")

pos0.95 <- results$pos[which(results$tdg09 >= 0.95)]
points_col_vec <- rep("#cccccc", length(pos0.95))

for (i in pos0.95) {
  if (i %in% conv_sub_sites$pos_cln) {
    points_col_vec[which(pos0.95 == i)] <- conv_sub_sites$col[which(conv_sub_sites$pos_cln == i)]
  }
}

points(which(results$tdg09 >= 0.95),
       results$tdg09[results$tdg09 >= 0.95],
       col="#cccccc",
       pch=20)

points(pos0.95,
       results[which(results$pos %in% pos0.95),]$tdg09,
       col=points_col_vec,
       pch=20)

abline(h=0.95,
       lty=3)

# topo

plot(results$topo, type="h", col=col_vec,
     frame.plot = F, xlim=c(1,274), ylim=c(0,1.0),
     ylab="", xlab="", xaxt="n")

pos0.95 <- results$pos[which(results$topo >= 0.95)]
points_col_vec <- rep("#cccccc", length(pos0.95))

for (i in pos0.95) {
  if (i %in% conv_sub_sites$pos_cln) {
    points_col_vec[which(pos0.95 == i)] <- conv_sub_sites$col[which(conv_sub_sites$pos_cln == i)]
  }
}

points(which(results$topo >= 0.95),
       results$topo[results$topo >= 0.95],
       col="#cccccc",
       pch=20)

points(pos0.95,
       results[which(results$pos %in% pos0.95),]$topo,
       col=points_col_vec,
       pch=20)

abline(h=0.95,
       lty=3)


# venn diagrams

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
