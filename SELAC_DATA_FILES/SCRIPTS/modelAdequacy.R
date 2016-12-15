
library(selac)

PlotAdequacyResults <- function(adequate.obj.gtr, adequate.obj.selac.nog, adequate.obj.selac.wg, alpha=.1, known.functionality=0.9293696, file.name="modelAdequacy.pdf"){
    prop.interval <- seq(0,1 , by=0.05)
    pdf(file.name)
    plot(prop.interval, adequate.obj.gtr[[1]], ylab="", xlab="", pch=19, xlim=c(0,1), ylim=c(.5, 1), axes=FALSE, col=0)

    for(rep.index in 1:100){
        lines(prop.interval, pp[[rep.index]], col=add.alpha("blue", alpha=alpha))
    }

    for(rep.index in 1:100){
        lines(prop.interval, adequate.obj.selac.nog[[rep.index]], col=add.alpha("red", alpha=alpha))
    }

    for(rep.index in 1:100){
        lines(prop.interval, adequate.obj.selac.wg[[rep.index]], col=add.alpha("green", alpha=alpha))
    }

    title(xlab = "Proportion of true edge length", line=2)
    title(ylab = "Functionality", line=3)
    par(tck=.01)
    axis(2, at = seq(.5, 1, by = .1), las =1, lwd=1, labels=TRUE, mgp=c(1,.5,0))
    axis(1, at = seq(0, 1, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(1,.5,0))
    abline(h=known.functionality, lty=2)
    dev.off()

}

## Adequacy of selac + Gamma ##
library(selac)
load("yeastSalRokSelacNoG.Rdata")
selac.nog <- result
load("yeastSalRokSelacGamma.Rdata")
selac.wg <- result
load("yeastSalRokSelacGTRG.Rdata")
gtr.g <- result


print("Doing Scer")
## Brewer's yeast adequacy
pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="selac", model.to.simulate.under="selac", selac.obj.to.reconstruct=selac.wg, selac.obj.to.simulate=selac.wg, aa.optim.input=NULL, taxon.to.drop=1, partition.number=53, numcode=1)
save(pp, file="adequacy_Scer_reconSelacWG_simSelacWG.Rsave")

pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="selac", model.to.simulate.under="gtr", selac.obj.to.reconstruct=selac.wg, selac.obj.to.simulate=gtr.g, aa.optim.input=NULL, taxon.to.drop=1, partition.number=53, numcode=1)
save(pp, file="adequacy_Scer_reconSelacWG_simGTRG.Rsave")

pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="gtr", model.to.simulate.under="selac", selac.obj.to.reconstruct=gtr.g, selac.obj.to.simulate=selac.wg, aa.optim.input=NULL, taxon.to.drop=1, partition.number=53, numcode=1)
save(pp, file="adequacy_Scer_reconSelacGTRG_simSelacWG.Rsave")

pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="gtr", model.to.simulate.under="gtr", selac.obj.to.reconstruct=gtr.g, selac.obj.to.simulate=gtr.g, aa.optim.input=NULL, taxon.to.drop=1, partition.number=53, numcode=1, for.gtr.only=selac.wg)
save(pp, file="adequacy_Scer_reconGTRG_simGTRG.Rsave")

###############################

print("Doing Scas")

## Scas adequacy
pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="selac", model.to.simulate.under="selac", selac.obj.to.reconstruct=selac.wg, selac.obj.to.simulate=selac.wg, aa.optim.input=NULL, taxon.to.drop=5, partition.number=53, numcode=1)
save(pp, file="adequacy_Scas_reconSelacWG_simSelacWG.Rsave")

pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="selac", model.to.simulate.under="gtr", selac.obj.to.reconstruct=selac.wg, selac.obj.to.simulate=gtr.g, aa.optim.input=NULL, taxon.to.drop=5, partition.number=53, numcode=1)
save(pp, file="adequacy_Scas_reconSelacWG_simGTRG.Rsave")

pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="gtr", model.to.simulate.under="selac", selac.obj.to.reconstruct=gtr.g, selac.obj.to.simulate=selac.wg, aa.optim.input=NULL, taxon.to.drop=5, partition.number=53, numcode=1)
save(pp, file="adequacy_Scas_reconSelacGTRG_simSelacWG.Rsave")

pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="gtr", model.to.simulate.under="gtr", selac.obj.to.reconstruct=gtr.g, selac.obj.to.simulate=gtr.g, aa.optim.input=NULL, taxon.to.drop=5, partition.number=53, numcode=1, for.gtr.only=selac.wg)
save(pp, file="adequacy_Scas_reconGTRG_simGTRG.Rsave")

###############################



###Get functionality for a single taxon:

load("yeastSalRokSelacGamma.Rdata")
selac.wg <- result




