
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

#The code necessary for evaulating functionality across ALL genes:
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
    functionality <- c(functionality, GetFunctionality(gene.length=length(aa.optim), aa.data=aa.data[taxon.to.do,-1], optimal.aa=aa.optim, alpha=result$mle.pars[gene.index,2], beta=result$mle.pars[gene.index,3], gamma=0.0003990333))
}
write.table(cbind(result$partitions,functionality), file="function_taxon6", quote=FALSE, sep="\t", row.names=FALSE)





######################################################################################################################################
######################################################################################################################################
### NEW YEAST EMPIRICAL SUMMARY
######################################################################################################################################
######################################################################################################################################
gg <- read.delim("/finalPhiEsts.tsv")

pdf("SelACwG_vs_Empirical_by_spp.pdf", width=8, height=8)
par(mfcol=c(2,2),mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))

plot(log(gg$Phi_gamma/gg$functionality_Spar), log(gg$Spar_RNA), axes=FALSE, xlab="", ylab="", ylim=c(1,6), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(1,6, by = 1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~phi[RNA-seq]), line=2.5)
title(xlab=expression(log~hat(phi)[SelAC]), line=2)
fit <- lm(log(gg$Spar_RNA)~log(gg$Phi_gamma/gg$functionality_Spar))
abline(fit)
actual.max <- 6-1
max.diff <- actual.max-6
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(R^2 == .(round(summary.stats$r.squared,2)))
text(-1, bottom.text, eq)
text(-1, 6, "S. paradoxus")
mtext("(a)",side=3, line=0, adj=0)

plot(log(gg$Phi_gamma/gg$functionality_Smik), log(gg$Smik_RNA), axes=FALSE, xlab="", ylab="", ylim=c(2,10), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(2,10, by = 2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~phi[RNA-seq]), line=2.5)
title(xlab=expression(log~hat(phi)[SelAC]), line=2)
fit <- lm(log(gg$Smik_RNA)~log(gg$Phi_gamma/gg$functionality_Smik))
abline(fit)
actual.max <- 10-2
max.diff <- actual.max-10
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(R^2 == .(round(summary.stats$r.squared,2)))
text(-1, bottom.text, eq)
text(-1, 10, "S. mikatae")
mtext("(c)",side=3, line=0, adj=0)


plot(log(gg$Phi_gamma/gg$functionality_Scer), log(gg$Scer_RNA), axes=FALSE, xlab="", ylab="", ylim=c(-2.4,1.8), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(-2.4,1.8, by = .6), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~phi[RNA-seq]), line=2.5)
title(xlab=expression(log~hat(phi)[SelAC]), line=2)
fit <- lm(log(gg$Scer_RNA)~log(gg$Phi_gamma/gg$functionality_Scer))
abline(fit)
actual.max <- 1.8--2.4
max.diff <- actual.max-1.8
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(R^2 == .(round(summary.stats$r.squared,2)))
text(-1, bottom.text, eq)
text(-1, 1.8, "S. cervisiae")
mtext("(b)",side=3, line=0, adj=0)


plot(log(gg$Phi_gamma/gg$functionality_Scas), log(gg$Scas_Microarray), axes=FALSE, xlab="", ylab="", ylim=c(-3.5,0.5), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(-3.5,0.5, by = 1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~phi[Microarray]), line=2.5)
title(xlab=expression(log~hat(phi)[SelAC]), line=2)
fit <- lm(log(gg$Scas_Microarray)~log(gg$Phi_gamma/gg$functionality_Scas))
abline(fit)
actual.max <- .5--3.5
max.diff <- actual.max-.5
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(R^2 == .(round(summary.stats$r.squared,2)))
text(-1, bottom.text, eq)
text(-1, .5, "S. castellii")
mtext("(d)",side=3, line=0, adj=0)

dev.off()




pdf("SelACwoutG_vs_Empirical_by_spp.pdf", width=8, height=8)
par(mfcol=c(2,2),mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))

plot(log(gg$Phi_nogamma/gg$functionality_Spar), log(gg$Spar_RNA), axes=FALSE, xlab="", ylab="", ylim=c(1,6), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(1,6, by = 1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~psi[RNA-seq]), line=2.5)
title(xlab=expression(log~hat(psi)[SelAC]), line=2)
fit <- lm(log(gg$Spar_RNA)~log(gg$Phi_nogamma/gg$functionality_Spar))
abline(fit)
actual.max <- 6-1
max.diff <- actual.max-6
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(R^2 == .(round(summary.stats$r.squared,2)))
text(-1, bottom.text, eq)
text(-1, 6, "S. paradoxus")
mtext("(a)",side=3, line=0, adj=0)

plot(log(gg$Phi_nogamma/gg$functionality_Smik), log(gg$Smik_RNA), axes=FALSE, xlab="", ylab="", ylim=c(2,10), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(2,10, by = 2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~psi[RNA-seq]), line=2.5)
title(xlab=expression(log~hat(psi)[SelAC]), line=2)
fit <- lm(log(gg$Smik_RNA)~log(gg$Phi_nogamma/gg$functionality_Smik))
abline(fit)
actual.max <- 10-2
max.diff <- actual.max-10
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(R^2 == .(round(summary.stats$r.squared,2)))
text(-1, bottom.text, eq)
text(-1, 10, "S. mikatae")
mtext("(c)",side=3, line=0, adj=0)


plot(log(gg$Phi_nogamma/gg$functionality_Scer), log(gg$Scer_RNA), axes=FALSE, xlab="", ylab="", ylim=c(-2.4,1.8), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(-2.4,1.8, by = .6), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~psi[RNA-seq]), line=2.5)
title(xlab=expression(log~hat(psi)[SelAC]), line=2)
fit <- lm(log(gg$Scer_RNA)~log(gg$Phi_nogamma/gg$functionality_Scer))
abline(fit)
actual.max <- 1.8--2.4
max.diff <- actual.max-1.8
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(R^2 == .(round(summary.stats$r.squared,2)))
text(-1, bottom.text, eq)
text(-1, 1.8, "S. cervisiae")
mtext("(b)",side=3, line=0, adj=0)


plot(log(gg$Phi_nogamma/gg$functionality_Scas), log(gg$Scas_Microarray), axes=FALSE, xlab="", ylab="", ylim=c(-3.5,0.5), xlim=c(-1.2,0), pch=19, cex=.75, main="")
par(tck=.01)
axis(2, at = seq(-3.5,0.5, by = 1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2,0, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab=expression(log~psi[Microarray]), line=2.5)
title(xlab=expression(log~hat(psi)[SelAC]), line=2)
fit <- lm(log(gg$Scas_Microarray)~log(gg$Phi_nogamma/gg$functionality_Scas))
abline(fit)
actual.max <- .5--3.5
max.diff <- actual.max-.5
top.text <- (actual.max*.875)-max.diff
bottom.text <- (actual.max*.825)-max.diff
text(-1, top.text, paste(round(fit$coefficients[2],3), "x + ", round(fit$coefficients[1],2), sep=""))
summary.stats <- summary(fit)
eq <- bquote(R^2 == .(round(summary.stats$r.squared,2)))
text(-1, bottom.text, eq)
text(-1, .5, "S. castellii")
mtext("(d)",side=3, line=0, adj=0)

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



