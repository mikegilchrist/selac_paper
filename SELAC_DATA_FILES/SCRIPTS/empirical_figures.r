
######################################################################################################################################
######################################################################################################################################
### CALCULATING FUNCTIONALITY OF A GENE
######################################################################################################################################
######################################################################################################################################

library(selac)

#That function we call to evaluate ONE gene:
FunctionalityCalculation <- function(site.pattern.counts, aa.data, optimal.aa, alpha, beta, gamma){
    aa.properties <- structure(c(0, 2.75, 1.38, 0.92, 0, 0.74, 0.58, 0, 0.33, 0, 0,
    1.33, 0.39, 0.89, 0.65, 1.42, 0.71, 0, 0.13, 0.2, 8.1, 5.5, 13,
    12.3, 5.2, 9, 10.4, 5.2, 11.3, 4.9, 5.7, 11.6, 8, 10.5, 10.5,
    9.2, 8.6, 5.9, 5.4, 6.2, 31, 55, 54, 83, 132, 3, 96, 111, 119,
    111, 105, 56, 32.5, 85, 124, 32, 61, 84, 170, 136), .Dim = c(20L, 3L), .Dimnames = list(c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"), c("c", "p", "v"))) #properties from Grantham paper
    aa.distances <- c()
    #Note using only the second row, because we are comparing empirical S. cervisae rates:
    for(site.index in 1:site.pattern.counts){
        aa.distances <- c(aa.distances, (1+((alpha*(aa.properties[aa.data[2,site.index],1] - aa.properties[optimal.aa[site.index],1])^2 + beta*(aa.properties[aa.data[2,site.index],2]-aa.properties[optimal.aa[site.index],2])^2+gamma*(aa.properties[aa.data[2,site.index],3]-aa.properties[optimal.aa[site.index],3])^2)^(1/2))))
    }
    functionality = 1/((1/site.pattern.counts) * sum(aa.distances))
    return(functionality)
}

pp <- selac:::GetSelacPhiCat(result, codon.data.path="")
tt <- selac:::LaguerreQuad(result$mle.pars[1,15], 4)

#The code necessary for evaulating functionality across ALL genes:
phy <- result$phy
functionality <- c()
taxon.to.do <- 1
for(gene.index in 1:100){
    yeast.gene <- read.dna(result$partitions[gene.index], format="fasta")
    yeast.gene <- as.list(as.matrix(cbind(yeast.gene)))
    chars <- selac:::DNAbinToCodonNumeric(yeast.gene)
    codon.data <- chars[phy$tip.label,]
    aa.data <- selac:::ConvertCodonNumericDataToAAData(codon.data, numcode=1)
    aa.optim <- result$aa.optim[[gene.index]]
    functionality <- c(functionality, GetFunctionality(gene.length=length(aa.optim), gp=tt[pp[[gene.index]]$indicator.by.site.weightedLik], aa.data=aa.data[taxon.to.do,-1], optimal.aa=aa.optim, alpha=result$mle.pars[gene.index,2], beta=result$mle.pars[gene.index,3], gamma=0.0003990333))
}
write.table(cbind(result$partitions,functionality), file="function_taxon1", quote=FALSE, sep="\t", row.names=FALSE)


phy <- result$phy
functionality <- c()
taxon.to.do <- 2
for(gene.index in 1:100){
    yeast.gene <- read.dna(result$partitions[gene.index], format="fasta")
    yeast.gene <- as.list(as.matrix(cbind(yeast.gene)))
    chars <- selac:::DNAbinToCodonNumeric(yeast.gene)
    codon.data <- chars[phy$tip.label,]
    aa.data <- selac:::ConvertCodonNumericDataToAAData(codon.data, numcode=1)
    aa.optim <- result$aa.optim[[gene.index]]
    functionality <- c(functionality, GetFunctionality(gene.length=length(aa.optim), gp=tt[pp[[gene.index]]$indicator.by.site.weightedLik], aa.data=aa.data[taxon.to.do,-1], optimal.aa=aa.optim, alpha=result$mle.pars[gene.index,2], beta=result$mle.pars[gene.index,3], gamma=0.0003990333))
}
write.table(cbind(result$partitions,functionality), file="function_taxon2", quote=FALSE, sep="\t", row.names=FALSE)


phy <- result$phy
functionality <- c()
taxon.to.do <- 3
for(gene.index in 1:100){
    yeast.gene <- read.dna(result$partitions[gene.index], format="fasta")
    yeast.gene <- as.list(as.matrix(cbind(yeast.gene)))
    chars <- selac:::DNAbinToCodonNumeric(yeast.gene)
    codon.data <- chars[phy$tip.label,]
    aa.data <- selac:::ConvertCodonNumericDataToAAData(codon.data, numcode=1)
    aa.optim <- result$aa.optim[[gene.index]]
    functionality <- c(functionality, GetFunctionality(gene.length=length(aa.optim), gp=tt[pp[[gene.index]]$indicator.by.site.weightedLik], aa.data=aa.data[taxon.to.do,-1], optimal.aa=aa.optim, alpha=result$mle.pars[gene.index,2], beta=result$mle.pars[gene.index,3], gamma=0.0003990333))
}
write.table(cbind(result$partitions,functionality), file="function_taxon3", quote=FALSE, sep="\t", row.names=FALSE)


phy <- result$phy
functionality <- c()
taxon.to.do <- 4
for(gene.index in 1:100){
    yeast.gene <- read.dna(result$partitions[gene.index], format="fasta")
    yeast.gene <- as.list(as.matrix(cbind(yeast.gene)))
    chars <- selac:::DNAbinToCodonNumeric(yeast.gene)
    codon.data <- chars[phy$tip.label,]
    aa.data <- selac:::ConvertCodonNumericDataToAAData(codon.data, numcode=1)
    aa.optim <- result$aa.optim[[gene.index]]
    functionality <- c(functionality, GetFunctionality(gene.length=length(aa.optim), gp=tt[pp[[gene.index]]$indicator.by.site.weightedLik], aa.data=aa.data[taxon.to.do,-1], optimal.aa=aa.optim, alpha=result$mle.pars[gene.index,2], beta=result$mle.pars[gene.index,3], gamma=0.0003990333))
}
write.table(cbind(result$partitions,functionality), file="function_taxon4", quote=FALSE, sep="\t", row.names=FALSE)


phy <- result$phy
functionality <- c()
taxon.to.do <- 5
for(gene.index in 1:100){
    yeast.gene <- read.dna(result$partitions[gene.index], format="fasta")
    yeast.gene <- as.list(as.matrix(cbind(yeast.gene)))
    chars <- selac:::DNAbinToCodonNumeric(yeast.gene)
    codon.data <- chars[phy$tip.label,]
    aa.data <- selac:::ConvertCodonNumericDataToAAData(codon.data, numcode=1)
    aa.optim <- result$aa.optim[[gene.index]]
    functionality <- c(functionality, GetFunctionality(gene.length=length(aa.optim), gp=tt[pp[[gene.index]]$indicator.by.site.weightedLik], aa.data=aa.data[taxon.to.do,-1], optimal.aa=aa.optim, alpha=result$mle.pars[gene.index,2], beta=result$mle.pars[gene.index,3], gamma=0.0003990333))
}
write.table(cbind(result$partitions,functionality), file="function_taxon5", quote=FALSE, sep="\t", row.names=FALSE)


phy <- result$phy
functionality <- c()
taxon.to.do <- 6
for(gene.index in 1:100){
    yeast.gene <- read.dna(result$partitions[gene.index], format="fasta")
    yeast.gene <- as.list(as.matrix(cbind(yeast.gene)))
    chars <- selac:::DNAbinToCodonNumeric(yeast.gene)
    codon.data <- chars[phy$tip.label,]
    aa.data <- selac:::ConvertCodonNumericDataToAAData(codon.data, numcode=1)
    aa.optim <- result$aa.optim[[gene.index]]
    functionality <- c(functionality, GetFunctionality(gene.length=length(aa.optim), gp=tt[pp[[gene.index]]$indicator.by.site.weightedLik], aa.data=aa.data[taxon.to.do,-1], optimal.aa=aa.optim, alpha=result$mle.pars[gene.index,2], beta=result$mle.pars[gene.index,3], gamma=0.0003990333))
}
write.table(cbind(result$partitions,functionality), file="function_taxon6", quote=FALSE, sep="\t", row.names=FALSE)


#gp=tt[pp[[gene.index]]$indicator.by.site.weightedLik]


######################################################################################################################################
######################################################################################################################################
### NEW YEAST EMPIRICAL SUMMARY
######################################################################################################################################
######################################################################################################################################
gg <- read.delim("finalPhiEstsNEW_includingROC.tsv")

pdf("SelACwG_vs_Empirical_by_spp.pdf", width=8, height=8)
par(mfcol=c(2,2),mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))

plot(log(gg$Psi_gamma/gg$functionality_Spar), log(gg$Spar_RNA), axes=FALSE, xlab="", ylab="", ylim=c(1,6), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(1,6, by = 1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~phi[RNA-seq]), line=2.5)
title(xlab=expression(log~hat(phi)[SelAC]), line=2)
fit <- lm(log(gg$Spar_RNA)~log(gg$Psi_gamma/gg$functionality_Spar))
abline(fit)
actual.max <- 6-1
max.diff <- actual.max-6
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(sqrt(summary.stats$r.squared),2)))
text(-1, bottom.text, eq)
text(-1, 6, expression(italic('S. paradoxus')))
mtext("(a)",side=3, line=0, adj=0)

plot(log(gg$Psi_gamma/gg$functionality_Smik), log(gg$Smik_RNA), axes=FALSE, xlab="", ylab="", ylim=c(2,10), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(2,10, by = 2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~phi[RNA-seq]), line=2.5)
title(xlab=expression(log~hat(phi)[SelAC]), line=2)
fit <- lm(log(gg$Smik_RNA)~log(gg$Psi_gamma/gg$functionality_Smik))
abline(fit)
actual.max <- 10-2
max.diff <- actual.max-10
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(sqrt(summary.stats$r.squared),2)))
text(-1, bottom.text, eq)
text(-1, 10, expression(italic('S. mikatae')))
mtext("(c)",side=3, line=0, adj=0)


plot(log(gg$Psi_gamma/gg$functionality_Scer), log(gg$Scer_RNA), axes=FALSE, xlab="", ylab="", ylim=c(-2.4,1.8), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(-2.4,1.8, by = .6), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~phi[RNA-seq]), line=2.5)
title(xlab=expression(log~hat(phi)[SelAC]), line=2)
fit <- lm(log(gg$Scer_RNA)~log(gg$Psi_gamma/gg$functionality_Scer))
abline(fit)
actual.max <- 1.8--2.4
max.diff <- actual.max-1.8
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(sqrt(summary.stats$r.squared),2)))
text(-1, bottom.text, eq)
text(-1, 1.8, expression(italic('S. cervisiae')))
mtext("(b)",side=3, line=0, adj=0)


plot(log(gg$Psi_gamma/gg$functionality_Scas), log(gg$Scas_Microarray), axes=FALSE, xlab="", ylab="", ylim=c(-3.5,0.5), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(-3.5,0.5, by = 1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~phi[Microarray]), line=2.5)
title(xlab=expression(log~hat(phi)[SelAC]), line=2)
fit <- lm(log(gg$Scas_Microarray)~log(gg$Psi_gamma/gg$functionality_Scas))
abline(fit)
actual.max <- .5--3.5
max.diff <- actual.max-.5
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(sqrt(summary.stats$r.squared),2)))
text(-1, bottom.text, eq)
text(-1, .5, expression(italic('S. castellii')))
mtext("(d)",side=3, line=0, adj=0)

dev.off()



gg <- read.delim("finalPhiEsts_including_ROC.tsv")

pdf("SelACwG_vs_ROC_by_spp.pdf", width=8, height=8)
par(mfcol=c(2,2),mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))

plot(log(gg$Psi_gamma/gg$functionality_Spar_wG), log(gg$ROC_Spar), axes=FALSE, xlab="", ylab="", ylim=c(-2,1.5), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(-2,1.5, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~phi[ROC]), line=2.5)
title(xlab=expression(log~hat(phi)[SelAC]), line=2)
fit <- lm(log(gg$ROC_Spar)~log(gg$Psi_gamma/gg$functionality_Spar))
abline(fit)
actual.max <- 1.5--2
max.diff <- actual.max-1.5
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(sqrt(summary.stats$r.squared),2)))
text(-1, bottom.text, eq)
text(-1, 1.5, expression(italic('S. paradoxus')))
mtext("(a)",side=3, line=0, adj=0)

plot(log(gg$Psi_gamma/gg$functionality_Smik), log(gg$ROC_Smik), axes=FALSE, xlab="", ylab="", ylim=c(-2,1.5), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(-2,1.5, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~phi[ROC]), line=2.5)
title(xlab=expression(log~hat(phi)[SelAC]), line=2)
fit <- lm(log(gg$ROC_Smik)~log(gg$Psi_gamma/gg$functionality_Smik))
abline(fit)
actual.max <- 1.5--2
max.diff <- actual.max-1.5
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(sqrt(summary.stats$r.squared),2)))
text(-1, bottom.text, eq)
text(-1, 1.5, expression(italic('S. mikatae')))
mtext("(c)",side=3, line=0, adj=0)


plot(log(gg$Psi_gamma/gg$functionality_Scer), log(gg$ROC_Scer), axes=FALSE, xlab="", ylab="", ylim=c(-2,1.5), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(-2,1.5, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~phi[ROC]), line=2.5)
title(xlab=expression(log~hat(phi)[SelAC]), line=2)
fit <- lm(log(gg$ROC_Scer)~log(gg$Psi_gamma/gg$functionality_Scer))
abline(fit)
actual.max <- 1.5--2
max.diff <- actual.max-1.5
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(sqrt(summary.stats$r.squared),2)))
text(-1, bottom.text, eq)
text(-1, 1.5, expression(italic('S. cerevisiae')))
mtext("(b)",side=3, line=0, adj=0)


plot(log(gg$Psi_gamma/gg$functionality_Scas), log(gg$ROC_Scas), axes=FALSE, xlab="", ylab="", ylim=c(-2,1.5), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(-2,1.5, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~phi[ROC]), line=2.5)
title(xlab=expression(log~hat(phi)[SelAC]), line=2)
fit <- lm(log(gg$ROC_Scas)~log(gg$Psi_gamma/gg$functionality_Scas))
abline(fit)
actual.max <- 1.5--2
max.diff <- actual.max-1.5
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(sqrt(summary.stats$r.squared),2)))
text(-1, bottom.text, eq)
text(-1, 1.5, expression(italic('S. castellii')))
mtext("(d)",side=3, line=0, adj=0)

dev.off()


### Check the inclusion of functionality:
fit_Spar_wFunc <- lm(log(gg$ROC_Spar)~log(gg$Psi_gamma/gg$functionality_Spar))
summary.stats <- summary(fit_Spar_wFunc)
round(sqrt(summary.stats$r.squared),2)
fit_Spar_noFunc <- lm(log(gg$ROC_Spar)~log(gg$Psi_gamma))
summary.stats <- summary(fit_Spar_noFunc)
round(sqrt(summary.stats$r.squared),2)

fit_Smik_wFunc <- lm(log(gg$ROC_Smik)~log(gg$Psi_gamma/gg$functionality_Smik))
summary.stats <- summary(fit_Smik_wFunc)
round(sqrt(summary.stats$r.squared),2)
fit_Smik_noFunc <- lm(log(gg$ROC_Smik)~log(gg$Psi_gamma))
summary.stats <- summary(fit_Smik_noFunc)
round(sqrt(summary.stats$r.squared),2)

fit_Scer_wFunc <- lm(log(gg$ROC_Scer)~log(gg$Psi_gamma/gg$functionality_Scer))
summary.stats <- summary(fit_Scer_wFunc)
round(sqrt(summary.stats$r.squared),2)
fit_Scer_noFunc <- lm(log(gg$ROC_Scer)~log(gg$Psi_gamma))
summary.stats <- summary(fit_Scer_noFunc)
round(sqrt(summary.stats$r.squared),2)

fit_Scas_wFunc <- lm(log(gg$ROC_Scas)~log(gg$Psi_gamma/gg$functionality_Scas))
summary.stats <- summary(fit_Scas_wFunc)
round(sqrt(summary.stats$r.squared),2)
fit_Scas_noFunc <- lm(log(gg$ROC_Scas)~log(gg$Psi_gamma))
summary.stats <- summary(fit_Scas_noFunc)
round(sqrt(summary.stats$r.squared),2)


gg <- read.delim("/finalPhiEsts.tsv")

pdf("Empirical_vs_ROC_by_spp.pdf", width=8, height=8)
par(mfcol=c(2,2),mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))

plot(log(gg$Spar_RNA), log(gg$ROC_Spar), axes=FALSE, xlab="", ylab="", xlim=c(1,6), ylim=c(-2,1.5), pch=19, cex=.75, main="")
par(tck=.01)
axis(1, at = seq(1,6, by = 1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(2, at = seq(-2,1.5, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(xlab=expression(log~phi[RNA-seq]), line=2.5)
title(ylab=expression(log~phi[ROC]), line=2)
fit <- lm(log(gg$ROC_Spar)~log(gg$Spar_RNA))
abline(fit)
actual.max <- 1.5--2
max.diff <- actual.max-1.5
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(2, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(sqrt(summary.stats$r.squared),2)))
text(2, bottom.text, eq)
text(2, 1.5, expression(italic('S. paradoxus')))
mtext("(a)",side=3, line=0, adj=0)

plot(log(gg$Smik_RNA), log(gg$ROC_Smik), axes=FALSE, xlab="", ylab="", xlim=c(2,10), ylim=c(-2,1.5), pch=19, cex=.75, main="")
par(tck=.01)
axis(1, at = seq(2,10, by = 2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(2, at = seq(-2,1.5, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(xlab=expression(log~phi[RNA-seq]), line=2.5)
title(ylab=expression(log~phi[ROC]), line=2)
fit <- lm(log(gg$ROC_Smik) ~ log(gg$Smik_RNA))
abline(fit)
actual.max <- 1.5--2
max.diff <- actual.max-1.5
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(4, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(sqrt(summary.stats$r.squared),2)))
text(4, bottom.text, eq)
text(4, 1.5, expression(italic('S. mikatae')))
mtext("(c)",side=3, line=0, adj=0)


plot(log(gg$Scer_RNA), log(gg$ROC_Scer), axes=FALSE, xlab="", ylab="", xlim=c(-2.4,1.8), ylim=c(-2,1.5), pch=19, cex=.75, main="")
par(tck=.01)
axis(1, at = seq(-2.4,1.8, by = .6), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(2, at = seq(-2,1.5, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(xlab=expression(log~phi[RNA-seq]), line=2.5)
title(ylab=expression(log~phi[ROC]), line=2)
fit <- lm(log(gg$ROC_Scer)~log(gg$Scer_RNA))
abline(fit)
actual.max <- 1.5--2
max.diff <- actual.max-1.5
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1.4, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(sqrt(summary.stats$r.squared),2)))
text(-1.4, bottom.text, eq)
text(-1.4, 1.5, expression(italic('S. cervisiae')))
mtext("(b)",side=3, line=0, adj=0)


plot(log(gg$Scas_Microarray), log(gg$ROC_Scas), axes=FALSE, xlab="", ylab="", xlim=c(-3.5,0.5), ylim=c(-2,1.5), pch=19, cex=.75, main="")
par(tck=.01)
axis(1, at = seq(-3.5,0.5, by = 1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(2, at = seq(-2,1.5, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(xlab=expression(log~phi[Microarray]), line=2.5)
title(ylab=expression(log~phi[ROC]), line=2)
fit <- lm(log(gg$ROC_Scas) ~ log(gg$Scas_Microarray))
abline(fit)
actual.max <- 1.5--2
max.diff <- actual.max-1.5
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-2.5, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(sqrt(summary.stats$r.squared),2)))
text(-2.5, bottom.text, eq)
text(-2.5, 1.5, expression(italic('S. castellii')))
mtext("(d)",side=3, line=0, adj=0)

dev.off()



gg <- read.delim("finalPhiEsts_including_ROC.tsv")

pdf("MutSelOmega_vs_Us_ROC_Scer_only.pdf", width=8, height=8)
par(mfcol=c(2,2),mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))

plot(log(gg$Omega), log(gg$Psi_gamma), axes=FALSE, xlab="", ylab="", xlim=c(-4.5,-1), ylim=c(-1.4,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(1, at = seq(-4.5,-1, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(2, at = seq(-1.4, 0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(xlab=expression(log~hat(omega)[FMutSel]), line=2.5)
title(ylab=expression(log~hat(psi)[SelAC]), line=2)
fit <- lm(log(gg$Psi_gamma)~log(gg$Omega))
abline(fit)
actual.max <- 0--1.2
max.diff <- actual.max-0
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1.5, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(-sqrt(summary.stats$r.squared),2)))
text(-1.5, bottom.text, eq)
mtext("(a)",side=3, line=0, adj=0)


plot(log(gg$Omega), log(gg$ROC_Scer), axes=FALSE, xlab="", ylab="", xlim=c(-4.5,-1), ylim=c(-2,1.5), pch=19, cex=.75, main="")
par(tck=.01)
axis(1, at = seq(-4.5,-1, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(2, at = seq(-2,1.5, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~phi[ROC]), line=2.5)
title(xlab=expression(log~hat(omega)[FMutSel]), line=2.5)
fit <- lm(log(gg$ROC_Scer) ~ log(gg$Omega))
abline(fit)
actual.max <- 1.5--2
max.diff <- actual.max-1.5
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1.5, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(-sqrt(summary.stats$r.squared),2)))
text(-1.5, bottom.text, eq)
mtext("(c)",side=3, line=0, adj=0)


plot(log(gg$Omega), log(gg$Scer_RNA), axes=FALSE, xlab="", ylab="", xlim=c(-4.5,-1), ylim=c(-2.4,1.8), pch=19, cex=.75, main="")
par(tck=.01)
axis(1, at = seq(-4.5,-1, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(2, at = seq(-2.4,1.8, by = .6), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~phi[RNA-seq]), line=2.5)
title(xlab=expression(log~hat(omega)[FMutSel]), line=2.5)
fit <- lm(log(gg$Scer_RNA)~log(gg$Omega))
abline(fit)
actual.max <- 1.8--2.4
max.diff <- actual.max-1.8
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1.4, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(r == .(round(-sqrt(summary.stats$r.squared),2)))
text(-1.4, bottom.text, eq)
mtext("(b)",side=3, line=0, adj=0)


dev.off()





#I forget why I have this code -- I think this was for a talk or something. Will keep for now.
pdf("fakeDiver.pdf")
plot(0:5, seq(0.5,3,by=.5),axes=FALSE, xlab="", ylab="", ylim=c(.5,3), xlim=c(0, 5), col=0, main="")
par(tck=.01)
axis(2, at = seq(.5,3, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,5, by = 1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab="ln(lineages)", line=2.5)
title(xlab="Age since clade origin", line=2)
dev.off()



