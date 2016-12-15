

######################################################################################################################################
######################################################################################################################################
### SIMULATION SUMMARY 5 GENES
######################################################################################################################################
######################################################################################################################################

library(viridis)
library(selac)

#setwd("/Users/jeremybeaulieu/Desktop/WORK_IN_PROGRESS/SELAC_FINAL_DATASETS/sims/nogammaSims/restarts5geneALL")

#Get best fit -- To deal with Newton wall time, I did individual restarts with different seeds.
for(rep.index in 1:100){
    load(paste("restart1/sim5genesRes.", rep.index, ".Rsave", sep=""))
    round1 <- result
    #load(paste("restart2/sim5genesRes2.", rep.index, ".Rsave", sep=""))
    #round2 <- result
    #load(paste("restart3/sim5genesRes3.", rep.index, ".Rsave", sep=""))
    #round3 <- result
    tmp.liks <- c(round1$loglik,-100000000, -100000000)
    best <- which.max(tmp.liks)
    if(best == 1){
        result <- round1
        save(result, file=paste("final5gene.", rep.index, ".Rsave", sep=""))
    }
    if(best == 2){
        result <- round2
        save(result, file=paste("final5gene.", rep.index, ".Rsave", sep=""))
    }
    if(best == 3){
        result <- round3
        save(result, file=paste("final5gene.", rep.index, ".Rsave", sep=""))
    }
}


#Now get results
res.phi <- c()
res.gran <- c()
res.freq <- c()
res.gtr <- c()
res.edge.lengths <- c()

for(rep.index in 1:100){
    load(paste("final5gene.", rep.index, ".Rsave", sep=""))
    res.phi <- rbind(res.phi, result$mle.pars[,1])
    res.gran <- rbind(res.gran, result$mle.pars[1,2:3])
    res.freqs <- rbind(res.freq, result$mle.pars[1,4:6])
    res.gtr <- rbind(res.gtr, result$mle.pars[1,4:14])
    res.edge.lengths <- rbind(res.edge.lengths, result$phy$edge.length)
}

obj<-NULL
obj$phi <- res.phi
obj$gran <- res.gran
obj$freq <- res.freq
obj$gtr <- res.gtr
obj$edge.lengths <- res.edge.lengths

load("../../yeast.practiceSELACunrest106.Rdata")
known.optimal.aa <- result$aa.optim
pars.mat <- result$mle.pars
rows.to.sample = order(pars.mat[,1])[seq(5,105,20)]
true.par.vec = NULL
true.par.vec$C.q.phi.Ne = c(pars.mat[rows.to.sample,1])
true.par.vec$gtr = c(pars.mat[1,4:14])

#Alpha, beta
pdf("Figure_5genes_alpha.beta_UNREST.pdf", width=8, height=8)
par(mfcol=c(2,2),mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))

hist(obj$gran[,1], axes=FALSE, xlab="", ylab="", xlim=c(.2,.7), ylim=c(0,30), main="", col="gray")
par(tck=.01)
axis(2, at = seq(0, 30, by = 5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(.2, .7, by = .1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab="Frequency", line=2.5)
title(xlab=expression(Estimated~g[alpha]), line=2)
abline(v=0.32339585, lwd=2, lty=2)

plot(NULL)

hist(obj$gran[,2], axes=FALSE, xlab="", ylab="", xlim=c(.05, .20), main="", ylim=c(0,30), col="gray")
par(tck=.01)
axis(1, at = seq(.05, .20, by = .025), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(2, at = seq(0, 30, by = 5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab="Frequency", line=2.5)
title(xlab=expression(Estimated~g[beta]), line=2)
abline(v=0.09785295, lwd=2, lty=2)

dev.off()


###Log Phi, GTR, Optimal AA

#Phi
known <- true.par.vec$C.q.phi.Ne
C=4
q=4e-7
Phi.Ne.q <- obj$phi / C
Phi.Ne <- Phi.Ne.q/q
Phi.est <- Phi.Ne/5e6

Phi.Ne.q <- known / C
Phi.Ne <- Phi.Ne.q/q
Phi.known <- Phi.Ne/5e6

pdf("Figure_5genes_All_UNREST.pdf")
par(mfcol=c(2,2),mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
color.scheme <- viridis(5)
plot(colMeans(log(Phi.est)), log(Phi.known[1:5]), axes=FALSE, xlab="", ylab="", xlim=c(-1.8,-.8), ylim=c(-1.8,-.8), pch=19, col=0)
abline(0,1, lty=2)
low.upp <- apply(log(Phi.est), 2, quantile, c(0.025,0.975))
#segments(colMeans(log(Phi.est)), log(Phi.known[1:5]), colMeans(log(Phi.est)), low.upp[2,])
#segments(colMeans(log(Phi.est)), log(Phi.known[1:5]), colMeans(log(Phi.est)), low.upp[1,])
segments(colMeans(log(Phi.est)), log(Phi.known[1:5]), low.upp[2,], log(Phi.known[1:5]))
segments(colMeans(log(Phi.est)), log(Phi.known[1:5]), low.upp[1,], log(Phi.known[1:5]))
points(colMeans(log(Phi.est)), log(Phi.known[1:5]), pch=19, col=color.scheme)
par(tck=.01)
axis(2, at = seq(-1.8,-.8, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.8,-.8, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(xlab=expression(log~hat(phi)), line=2.5)
title(ylab=expression(log~phi[true]), line=2)
mtext("(a)",side=3, line=0, adj=0)
#legend("topleft", c("Gene 1", "Gene 2", "Gene 3", "Gene 4", "Gene 5"), fill=color.scheme, box.col=0)
############################################################


## Optimal AA #
optimal.matches <- as.list(1:5)
load("../../yeast.practiceSELACunrest106.Rdata")
known.optimal.aa <- result$aa.optim
pars.mat <- result$mle.pars
rows.to.sample = order(pars.mat[,1])[seq(5,105,20)]
for(gene.index in 1:5){
    tmp <- c()
    for(rep.index in 1:100){
        load(paste("final5gene.", rep.index, ".Rsave", sep=""))
        tmp <- c(tmp, sum(result$aa.optim[[gene.index]] == known.optimal.aa[[rows.to.sample[gene.index]]])/length(result$aa.optim[[gene.index]]))
    }
    optimal.matches[[gene.index]] = tmp
}

means <- c()
low <- c()
upp <- c()
for(gene.index in 1:5){
    means <- c(means, mean(optimal.matches[[gene.index]]))
    tmp <- quantile(optimal.matches[[gene.index]], prob=c(0.025,0.975))
    low <- c(low, tmp[1])
    upp <- c(upp, tmp[2])
}
plot(1:5, means, axes=FALSE, xlab="", ylab="", ylim=c(0,1), xlim=c(0,6), col=viridis(5), type="h", lwd=30, lend=1)
par(tck=.01)
axis(2, at = seq(0,1, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,6, by = 1), las=1, lwd=1, labels=c("","1","2","3","4","5",""), mgp=c(.75,.5,0))
title(ylab="Proportion correct amino acid", line=2.5)
title(xlab="Gene index", line=2)
segments(1:5, means, 1:5, low)
segments(1:5, means, 1:5, upp)
mtext("(c)",side=3, line=0, adj=0)
############################################################


#GTR params log
library(viridis)
gtr.known <- true.par.vec$gtr
gtr.est <- obj$gtr

color.scheme <- viridis(6)[1]
plot(colMeans(log(gtr.est)), log(gtr.known), axes=FALSE, xlab="", ylab="", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5),pch=19, col=0)
low.upp <- apply(log(gtr.est), 2, quantile, c(0.025,0.975))
segments(colMeans(log(gtr.est)), log(gtr.known), low.upp[2,], log(gtr.known))
segments(colMeans(log(gtr.est)), log(gtr.known), low.upp[1,], log(gtr.known))
points(colMeans(log(gtr.est)), log(gtr.known), pch=19,  col=color.scheme)

par(tck=.01)
axis(2, at = seq(-1.5,1.5, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.5,1.5, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab="log True GTR parameters", line=2.5)
title(xlab="log Estimated GTR parameters", line=2)
abline(0,1, lty=2)
mtext("(b)",side=3, line=0, adj=0)
#legend("topleft", c("Gene 1", "Gene 2", "Gene 3", "Gene 4", "Gene 5"), fill=color.scheme, box.col=0)
#dev.off()

load("../../yeast.practiceSELACunrest106.Rdata")
tree <- result$phy
pp <- colMeans(obj$edge.lengths)
color.scheme <- viridis(6)[1]
plot(log(pp), log(tree$edge.length), axes=FALSE, xlab="", ylab="", xlim=c(-1,2), ylim=c(-1,2), pch=19, col=0)
par(tck=.01)
axis(2, at = seq(-1,2, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1,2, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab="log True edge lengths", line=2.5)
title(xlab="log Estimated edge lengths", line=2)
low.upp <- apply(log(obj$edge.lengths), 2, quantile, c(0.025,0.975))
segments(log(pp), log(tree$edge.length), log(pp), low.upp[2,])
segments(log(pp), log(tree$edge.length), log(pp), low.upp[1,])
points(log(pp), log(tree$edge.length), pch=19,  col=color.scheme)
abline(0,1, lty=2)
mtext("(d)",side=3, line=0, adj=0)

dev.off()


#Tree
pdf("Figure_5genes_EdgeLengths_UNREST.pdf")
par(mfrow=(c(1,2)))
load("../../yeast.practiceSELACunrest106.Rdata")
plot(result$phy, main="True tree")
for(rep.index in 1:100){
    load(paste("final5gene.", rep.index, ".Rsave", sep=""))
    if(rep.index == 100){
        plot(result$phy, edge.color=rgb(0, 0, 0, .1), show.tip.label=FALSE, main="Estimates")
    }else{
        plot(result$phy, edge.color=rgb(0, 0, 0, .1), show.tip.label=FALSE)
    }
    par(new=TRUE)
}
dev.off()



######################################################################################################################################
######################################################################################################################################
### SIMULATION SUMMARY 5 GENES with gamma
######################################################################################################################################
######################################################################################################################################


#Get best fit -- To deal with Newton wall time, I did individual restarts with different seeds.
files <- system(paste("ls -1 ", "*.Rsave", sep=""), intern=TRUE)
count <- 1
for(rep.index in seq(1,250,by=5)){
    load(files[rep.index])
    round1 <- result

    load(files[rep.index+1])
    round2 <- result

    load(files[rep.index+2])
    round3 <- result

    load(files[rep.index+3])
    round4 <- result

    load(files[rep.index+4])
    round5 <- result
    
    tmp.liks <- c(round1$loglik, round2$loglik, round3$loglik, round4$loglik, round5$loglik)
    best <- which.max(tmp.liks)
    if(best == 1){
        result <- round1
        save(result, file=paste("final5gene.", count, ".Rsave", sep=""))
    }
    if(best == 2){
        result <- round2
        save(result, file=paste("final5gene.", count, ".Rsave", sep=""))
    }
    if(best == 3){
        result <- round3
        save(result, file=paste("final5gene.", count, ".Rsave", sep=""))
    }
    if(best == 4){
        result <- round4
        save(result, file=paste("final5gene.", count, ".Rsave", sep=""))
    }
    if(best == 5){
        result <- round5
        save(result, file=paste("final5gene.", count, ".Rsave", sep=""))
    }
    count <- count + 1
}


#Now get results
res.phi <- c()
res.gran <- c()
res.freq <- c()
res.gtr <- c()
res.gammashape <- c()
res.edge.lengths <- c()

for(rep.index in 1:50){
    load(paste("final5gene.", rep.index, ".Rsave", sep=""))
    res.phi <- rbind(res.phi, result$mle.pars[,1])
    res.gran <- rbind(res.gran, result$mle.pars[1,2:3])
    res.freqs <- rbind(res.freq, result$mle.pars[1,4:6])
    res.gtr <- rbind(res.gtr, result$mle.pars[1,4:14])
    res.gammashape <- rbind(res.gammashape, result$mle.pars[1,15])
    res.edge.lengths <- rbind(res.edge.lengths, result$phy$edge.length)
}

obj<-NULL
obj$phi <- res.phi
obj$gran <- res.gran
obj$freq <- res.freq
obj$gtr <- res.gtr
obj$gammashape <- res.gammashape
obj$edge.lengths <- res.edge.lengths

load("../../yeast.practiceSELACunrestgamma106.Rdata")
known.optimal.aa <- result$aa.optim
pars.mat <- result$mle.pars
rows.to.sample = order(pars.mat[,1])[seq(5,105,20)]
true.par.vec = NULL
true.par.vec$C.q.phi.Ne = c(pars.mat[rows.to.sample,1])
true.par.vec$gtr = c(pars.mat[1,4:14])
true.par.vec$gammashape = c(pars.mat[1,15])

#Alpha, beta
pdf("Figure_5genes_alpha.beta_gammashape_UNREST_WITHGAMMA.pdf", width=8, height=8)
par(mfcol=c(2,2),mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))

hist(obj$gran[,1], axes=FALSE, xlab="", ylab="", xlim=c(.2,.7), ylim=c(0,15), breaks=20, main="", col="gray")
par(tck=.01)
axis(2, at = seq(0, 15, by = 2.5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(.2, .7, by = .1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab="Frequency", line=2.5)
title(xlab=expression(Estimated~alpha[c]), line=2)
abline(v=0.4013160, lwd=2, lty=2)

hist(obj$gammashape, axes=FALSE, xlab="", ylab="", xlim=c(0, 2.5), main="", ylim=c(0,15), breaks=20, col="gray")
par(tck=.01)
axis(1, at = seq(0, 2.5, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(2, at = seq(0, 15, by = 2.5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab="Frequency", line=2.5)
title(xlab=expression(Estimated~alpha[g]), line=2)
abline(v=1.98, lwd=2, lty=2)

hist(obj$gran[,2], axes=FALSE, xlab="", ylab="", xlim=c(.05, .20), main="", ylim=c(0,15), breaks=20, col="gray")
par(tck=.01)
axis(1, at = seq(.05, .20, by = .025), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(2, at = seq(0, 15, by = 2.5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab="Frequency", line=2.5)
title(xlab=expression(Estimated~alpha[p]), line=2)
abline(v=0.1063492, lwd=2, lty=2)

dev.off()


###Log Phi, GTR, Optimal AA

#Phi
known <- true.par.vec$C.q.phi.Ne
C=4
q=4e-7
Phi.Ne.q <- obj$phi / C
Phi.Ne <- Phi.Ne.q/q
Phi.est <- Phi.Ne/5e6

Phi.Ne.q <- known / C
Phi.Ne <- Phi.Ne.q/q
Phi.known <- Phi.Ne/5e6

pdf("Figure_5genes_All_UNREST_WITHGAMMA.pdf")
par(mfcol=c(2,2),mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
color.scheme <- viridis(5)
plot(colMeans(log(Phi.est)), log(Phi.known[1:5]), axes=FALSE, xlab="", ylab="", xlim=c(-1.2,-.2), ylim=c(-1.2,-.2), pch=19, col=0)
abline(0,1, lty=2)
low.upp <- apply(log(Phi.est), 2, quantile, c(0.025,0.975))
#segments(colMeans(log(Phi.est)), log(Phi.known[1:5]), colMeans(log(Phi.est)), low.upp[2,])
#segments(colMeans(log(Phi.est)), log(Phi.known[1:5]), colMeans(log(Phi.est)), low.upp[1,])
segments(colMeans(log(Phi.est)), log(Phi.known[1:5]), low.upp[2,], log(Phi.known[1:5]))
segments(colMeans(log(Phi.est)), log(Phi.known[1:5]), low.upp[1,], log(Phi.known[1:5]))
points(colMeans(log(Phi.est)), log(Phi.known[1:5]), pch=19, col=color.scheme)
par(tck=.01)
axis(2, at = seq(-1.2, -.2, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.2, -.2, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(xlab=expression(log~hat(phi)), line=2.5)
title(ylab=expression(log~phi[true]), line=2)
mtext("(a)",side=3, line=0, adj=0)
#legend("topleft", c("Gene 1", "Gene 2", "Gene 3", "Gene 4", "Gene 5"), fill=color.scheme, box.col=0)
############################################################


## Optimal AA #
optimal.matches <- as.list(1:5)
load("../../yeast.practiceSELACunrestgamma106.Rdata")
known.optimal.aa <- result$aa.optim
pars.mat <- result$mle.pars
rows.to.sample = order(pars.mat[,1])[seq(5,105,20)]
for(gene.index in 1:5){
    tmp <- c()
    for(rep.index in 1:50){
        load(paste("final5gene.", rep.index, ".Rsave", sep=""))
        tmp <- c(tmp, sum(result$aa.optim[[gene.index]] == known.optimal.aa[[rows.to.sample[gene.index]]])/length(result$aa.optim[[gene.index]]))
    }
    optimal.matches[[gene.index]] = tmp
}

means <- c()
low <- c()
upp <- c()
for(gene.index in 1:5){
    means <- c(means, mean(optimal.matches[[gene.index]]))
    tmp <- quantile(optimal.matches[[gene.index]], prob=c(0.025,0.975))
    low <- c(low, tmp[1])
    upp <- c(upp, tmp[2])
}
plot(1:5, means, axes=FALSE, xlab="", ylab="", ylim=c(0,1), xlim=c(0,6), col=viridis(5), type="h", lwd=30, lend=1)
par(tck=.01)
axis(2, at = seq(0,1, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,6, by = 1), las=1, lwd=1, labels=c("","1","2","3","4","5",""), mgp=c(.75,.5,0))
title(ylab="Proportion correct amino acid", line=2.5)
title(xlab="Gene index", line=2)
segments(1:5, means, 1:5, low)
segments(1:5, means, 1:5, upp)
mtext("(c)",side=3, line=0, adj=0)
############################################################


#GTR params log
library(viridis)
gtr.known <- true.par.vec$gtr
gtr.est <- obj$gtr

color.scheme <- viridis(6)[1]
plot(colMeans(log(gtr.est)), log(gtr.known), axes=FALSE, xlab="", ylab="", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5),pch=19, col=0)
low.upp <- apply(log(gtr.est), 2, quantile, c(0.025,0.975))
segments(colMeans(log(gtr.est)), log(gtr.known), low.upp[2,], log(gtr.known))
segments(colMeans(log(gtr.est)), log(gtr.known), low.upp[1,], log(gtr.known))
points(colMeans(log(gtr.est)), log(gtr.known), pch=19,  col=color.scheme)

par(tck=.01)
axis(2, at = seq(-1.5,1.5, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1.5,1.5, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab="log True GTR parameters", line=2.5)
title(xlab="log Estimated GTR parameters", line=2)
abline(0,1, lty=2)
mtext("(b)",side=3, line=0, adj=0)
#legend("topleft", c("Gene 1", "Gene 2", "Gene 3", "Gene 4", "Gene 5"), fill=color.scheme, box.col=0)
#dev.off()

load("../../yeast.practiceSELACunrestgamma106.Rdata")
tree <- result$phy
pp <- colMeans(obj$edge.lengths)
color.scheme <- viridis(6)[1]
plot(log(pp), log(tree$edge.length), axes=FALSE, xlab="", ylab="", xlim=c(-1,2), ylim=c(-1,2), pch=19, col=0)
par(tck=.01)
axis(2, at = seq(-1,2, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-1,2, by = .5), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(ylab="log True edge lengths", line=2.5)
title(xlab="log Estimated edge lengths", line=2)
low.upp <- apply(log(obj$edge.lengths), 2, quantile, c(0.025,0.975))
segments(log(pp), log(tree$edge.length), log(pp), low.upp[2,])
segments(log(pp), log(tree$edge.length), log(pp), low.upp[1,])
points(log(pp), log(tree$edge.length), pch=19,  col=color.scheme)
abline(0,1, lty=2)
mtext("(d)",side=3, line=0, adj=0)

dev.off()


#Tree
pdf("Figure_5genes_EdgeLengths_UNREST.pdf")
par(mfrow=(c(1,2)))
load("../../yeast.practiceSELACunrest106.Rdata")
plot(result$phy, main="True tree")
for(rep.index in 1:100){
    load(paste("final5gene.", rep.index, ".Rsave", sep=""))
    if(rep.index == 100){
        plot(result$phy, edge.color=rgb(0, 0, 0, .1), show.tip.label=FALSE, main="Estimates")
    }else{
        plot(result$phy, edge.color=rgb(0, 0, 0, .1), show.tip.label=FALSE)
    }
    par(new=TRUE)
}
dev.off()

