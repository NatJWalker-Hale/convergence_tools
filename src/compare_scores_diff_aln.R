setwd("~/Dropbox/cary_projects/DODA/new_conv_2021/relaxed/20220117_restart/recon_topo/mafft/br_part_1st_origin/")

siteCorres <- read.table("site_corres_to_strict_recon_topo_prank_br_len_cln_aln.tsv", header=TRUE)

colnames(siteCorres) <- c("master", "comp")

# remove sites that are unaligned to either alignment

siteCorres <- siteCorres[which(siteCorres$master != 0),]
siteCorres <- siteCorres[which(siteCorres$comp != 0),]

masterScore <- read.csv("~/Dropbox/cary_projects/DODA/new_conv_2021/strict/recon_topo/prank/brlen_cln_aln/br_part_1st_origin/all_output.csv", header=TRUE)[,2:7]
compScore <- read.csv("all_output.csv", header=TRUE)

masterScoreIn <- masterScore[which(masterScore$pos %in% siteCorres$master),2:6]
compScoreIn <- compScore[which(compScore$pos %in% siteCorres$comp),2:6]

colnames(masterScoreIn) <- gsub("$", "_1", colnames(masterScoreIn))
colnames(compScoreIn) <- gsub("$", "_2", colnames(compScoreIn))

bothScores <- cbind(masterScoreIn, compScoreIn)

bothScores$diffsel_1_rank <- rank(bothScores$diffsel_1, ties.method = "average")
bothScores$diffsel_2_rank <- rank(bothScores$diffsel_2, ties.method = "average")

bothScores$msd_1_rank <- rank(bothScores$msd_1, ties.method = "average")
bothScores$msd_2_rank <- rank(bothScores$msd_2, ties.method = "average")

bothScores$pcoc_1_rank <- rank


plot(bothScores$diffsel_2 ~ bothScores$diffsel_1)
plot(bothScores$diffsel_2_rank ~ bothScores$diffsel_1_rank)
cor(bothScores$diffsel_1, bothScores$diffsel_2, method="kendall")

plot(bothScores$msd_2 ~ bothScores$msd_1)
plot(bothScores$msd_2_rank ~ bothScores$msd_1_rank)
cor(bothScores$msd_1, bothScores$msd_2, method="kendall")
