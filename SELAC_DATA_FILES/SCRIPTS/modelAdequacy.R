
library(selac)
library(viridis)

PlotAdequacyResults <- function(adequate.obj.gtr, adequate.obj.mutsel, adequate.obj.selac.wg, alpha=.2, known.functionality=0.88, file.name="modelAdequacyScer.pdf", taxon.names="S. cerevisiae"){
    prop.interval <- seq(0,1 , by=0.05)
    #pdf(file.name)
    cols <- viridis(3)
    plot(prop.interval, adequate.obj.gtr[[1]], ylab="", xlab="", pch=19, xlim=c(0,1), ylim=c(0, 1), axes=FALSE, col=0)

    for(rep.index in 1:100){
        lines(prop.interval, pp[[rep.index]], col=selac:::add.alpha(cols[1], alpha=alpha))
    }

    for(rep.index in 1:100){
        lines(prop.interval, adequate.obj.mutsel[[rep.index]], col=selac:::add.alpha(cols[2], alpha=alpha))
    }

    for(rep.index in 1:100){
        lines(prop.interval, adequate.obj.selac.wg[[rep.index]], col=selac:::add.alpha(cols[3], alpha=alpha))
    }

    title(xlab = "Proportion of true edge length", line=2)
    title(ylab = "Functionality", line=3)
    par(tck=.01)
    axis(2, at = seq(0, 1, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(1,.5,0))
    axis(1, at = seq(0, 1, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(1,.5,0))
    abline(h=known.functionality, lty=2)
    text(.15, .4, "SelAC", col=cols[3])
    text(.15, .35, "FMutSel", col=cols[2])
    text(.15, .3, expression(GTR+Gamma), col=cols[1])
    text(.15, 1, taxon.names)
    #dev.off()

}

pdf("modelAdequacyBOTH.pdf", width=8, height=8)
par(mfcol=c(2,2),mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
plot(NULL)
plot(NULL)
load("adequacy_Scer_reconSelacWG_simSelacWG.Rsave")
selac.wg <- pp
load("adequacy_Scer_reconSelacWG_simFMutSel.Rsave")
selac.mutsel <- pp
load("adequacy_Scer_reconSelacWG_simGTRG.Rsave")
selac.gtr <- pp
PlotAdequacyResults(adequate.obj.gtr=selac.gtr, adequate.obj.mutsel=selac.mutsel, adequate.obj.selac.wg=selac.wg, known.functionality=0.957376007, file.name="modelAdequacyScer.pdf",taxon.names="S. cerevisiae")
load("adequacy_Scas_reconSelacWG_simSelacWG.Rsave")
selac.wg <- pp
load("adequacy_Scas_reconSelacWG_simFMutSel.Rsave")
selac.mutsel <- pp
load("adequacy_Scas_reconSelacWG_simGTRG.Rsave")
selac.gtr <- pp
PlotAdequacyResults(adequate.obj.gtr=selac.gtr, adequate.obj.mutsel=selac.mutsel, adequate.obj.selac.wg=selac.wg, known.functionality=0.880932978, file.name="modelAdequacyScas.pdf",taxon.names="S. castellii")
dev.off()



## Adequacy of selac + Gamma ##
library(selac)
load("yeastSalRokSelacUNRESTnoG_NEW_.Rdata")
selac.nog <- result
load("yeastSalRokSelacUNRESTgamma_NEW_.Rdata")
selac.wg <- result
load("yeastSalRokGTRG.Rdata")
gtr.g <- result
load("yeastSalRokSelacFMutSel.Rdata")
fmutsel <- result
setwd("../DATA")

print("Doing Scer")
## Brewer's yeast adequacy
pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="selac", model.to.simulate.under="selac", selac.obj.to.reconstruct=selac.wg, selac.obj.to.simulate=selac.wg, aa.optim.input=NULL, taxon.to.drop=1, partition.number=53, numcode=1)
save(pp, file="adequacy_Scer_reconSelacWG_simSelacWG.Rsave")

pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="selac", model.to.simulate.under="gtr", selac.obj.to.reconstruct=selac.wg, selac.obj.to.simulate=gtr.g, aa.optim.input=NULL, taxon.to.drop=1, partition.number=53, numcode=1)
save(pp, file="adequacy_Scer_reconSelacWG_simGTRG.Rsave")

#pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="gtr", model.to.simulate.under="selac", selac.obj.to.reconstruct=gtr.g, selac.obj.to.simulate=selac.wg, aa.optim.input=NULL, taxon.to.drop=1, partition.number=53, numcode=1)
#save(pp, file="adequacy_Scer_reconSelacGTRG_simSelacWG.Rsave")

#pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="gtr", model.to.simulate.under="gtr", selac.obj.to.reconstruct=gtr.g, selac.obj.to.simulate=gtr.g, aa.optim.input=NULL, taxon.to.drop=1, partition.number=53, numcode=1, for.gtr.only=selac.wg)
#save(pp, file="adequacy_Scer_reconGTRG_simGTRG.Rsave")

pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="selac", model.to.simulate.under="fmutsel", selac.obj.to.reconstruct=selac.wg, selac.obj.to.simulate=fmutsel, aa.optim.input=NULL, taxon.to.drop=1, partition.number=53, numcode=1)
save(pp, file="adequacy_Scer_reconSelacWG_simFMutSel.Rsave")


###############################

print("Doing Scas")

## Scas adequacy
pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="selac", model.to.simulate.under="selac", selac.obj.to.reconstruct=selac.wg, selac.obj.to.simulate=selac.wg, aa.optim.input=NULL, taxon.to.drop=5, partition.number=53, numcode=1)
save(pp, file="adequacy_Scas_reconSelacWG_simSelacWG.Rsave")

pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="selac", model.to.simulate.under="gtr", selac.obj.to.reconstruct=selac.wg, selac.obj.to.simulate=gtr.g, aa.optim.input=NULL, taxon.to.drop=5, partition.number=53, numcode=1)
save(pp, file="adequacy_Scas_reconSelacWG_simGTRG.Rsave")

#pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="gtr", model.to.simulate.under="selac", selac.obj.to.reconstruct=gtr.g, selac.obj.to.simulate=selac.wg, aa.optim.input=NULL, taxon.to.drop=5, partition.number=53, numcode=1)
#save(pp, file="adequacy_Scas_reconSelacGTRG_simSelacWG.Rsave")

#pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="gtr", model.to.simulate.under="gtr", selac.obj.to.reconstruct=gtr.g, selac.obj.to.simulate=gtr.g, aa.optim.input=NULL, taxon.to.drop=5, partition.number=53, numcode=1, for.gtr.only=selac.wg)
#save(pp, file="adequacy_Scas_reconGTRG_simGTRG.Rsave")

pp <- GetAdequateManyReps(nreps=100, n.cores=4, model.to.reconstruct.under="selac", model.to.simulate.under="fmutsel", selac.obj.to.reconstruct=selac.wg, selac.obj.to.simulate=fmutsel, aa.optim.input=NULL, taxon.to.drop=5, partition.number=53, numcode=1)
save(pp, file="adequacy_Scas_reconSelacWG_simFMutSel.Rsave")


###############################





