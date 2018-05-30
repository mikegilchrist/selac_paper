library(lmodel2)

selac.exp.dat <- read.table(file = "finalPhiEsts.tsv", header=T, as.is=T, sep="\t")
######################
# GET PHI VALUE DATAFRAME
######################
selac.exp.phi <- data.frame(Scer_gamma_phi=selac.exp.dat$Psi_gamma/selac.exp.dat$functionality_Scer_wG,
                            Scer_nogamma_phi=selac.exp.dat$Psi_nogamma/selac.exp.dat$functionality_Scer_wG, 
                            Scas_gamma_phi=selac.exp.dat$Psi_gamma/selac.exp.dat$functionality_Scas_wG,
                            Scas_nogamma_phi=selac.exp.dat$Psi_nogamma/selac.exp.dat$functionality_Scas_wG, 
                            Smik_gamma_phi=selac.exp.dat$Psi_gamma/selac.exp.dat$functionality_Smik_wG,
                            Smik_nogamma_phi=selac.exp.dat$Psi_nogamma/selac.exp.dat$functionality_Smik_wG,
                            Skud_gamma_phi=selac.exp.dat$Psi_gamma/selac.exp.dat$functionality_Skud_wG,
                            Skud_nogamma_phi=selac.exp.dat$Psi_nogamma/selac.exp.dat$functionality_Skud_wG,
                            Spar_gamma_phi=selac.exp.dat$Psi_gamma/selac.exp.dat$functionality_Spar_wG,
                            Spar_nogamma_phi=selac.exp.dat$Psi_nogamma/selac.exp.dat$functionality_Spar_wG,
                            Cgla_gamma_phi=selac.exp.dat$Psi_gamma/selac.exp.dat$functionality_Cgla_wG,
                            Cgla_nogamma_phi=selac.exp.dat$Psi_nogamma/selac.exp.dat$functionality_Cgla_wG)




######################
# COMAPRE PSI ESTIMATES, GAMMA v. NO-GAMMA
######################

pdf("selac_psi_gamma_v_nogamma.pdf", height=5, width=4.9)
par(mar=c(5,5,4,2)+0.1)
plot(selac.exp.dat$Psi_gamma, selac.exp.dat$Psi_nogamma, 
     xlab=expression(widehat(psi)[gamma]), ylab=expression(widehat(psi)),
     main = expression("Effects of +"~Gamma~" on SELAC"~widehat(psi)))
sma.model <- lmodel2(Psi_nogamma~Psi_gamma, data = selac.exp.dat, range.y=NULL, range.x=NULL, nperm=0)
line.res <- as.numeric(sma.model$regression.results[3,2:3])
abline(a = line.res[1], b = line.res[2])
legend("bottomright", bty="n", legend = substitute(paste(rho, "=", r), list(r=round(sma.model$r, digits = 3))))
dev.off()


######################
# COMAPRE PHI ESTIMATES, GAMMA v. NO-GAMMA
######################

plotPhiComp <- function(data, organism, full.name){
  gamma.org <- paste(organism, "_gamma_phi", sep="")  
  nogamma.org <- paste(organism, "_nogamma_phi", sep="")
  plot(data[[gamma.org]], data[[nogamma.org]],
       xlab=expression(widehat(phi)[Gamma]), ylab=expression(widehat(phi)),
       main = expression("Effects of +"~Gamma~" on SELAC"~widehat(phi)))
  formula.reg <- paste(organism, "_nogamma_phi~", organism, "_gamma_phi", sep="")
  sma.model <- lmodel2(formula.reg, data = data, range.y=NULL, range.x=NULL, nperm=0)
  line.res <- as.numeric(sma.model$regression.results[3,2:3])
  abline(a = line.res[1], b = line.res[2])
  legend("bottomright", bty="n", legend = substitute(paste(rho, "=", r), list(r=round(sma.model$r, digits = 3))))
  legend("topleft", bty="n", legend = full.name, text.font=3)
}

pdf("selac_phi_gamma_v_nogamma.pdf", height=10, width=7)
par(mfrow=c(3,2), mar=c(4,5,3,1)+0.1)
plotPhiComp(selac.exp.phi, "Scer", "S. cereviciae")
plotPhiComp(selac.exp.phi, "Scas", "S. castellii")
plotPhiComp(selac.exp.phi, "Skud", "S. kudriavzevii")
plotPhiComp(selac.exp.phi, "Smik", "S. mikatae")
plotPhiComp(selac.exp.phi, "Spar", "S. paradoxus")
plotPhiComp(selac.exp.phi, "Cgla", "C. glabrata")
dev.off()

######################
# COMAPRE PHI ESTIMATES ACROSS SPECIES
######################

plotPhiCompSpecies <- function(data, organism1, organism2, full.name1, full.name2, gamma){
  
  org1.idx <- ifelse(gamma, 
                     paste(organism1, "_gamma_phi", sep=""),
                     paste(organism1, "_nogamma_phi", sep=""))
  org2.idx <- ifelse(gamma, 
                     paste(organism2, "_gamma_phi", sep=""),
                     paste(organism2, "_nogamma_phi", sep=""))
  
  xlabel <- substitute(paste(org, " ", widehat(phi)), list(org=full.name1))
  if(gamma){xlabel <- substitute(paste(org, " ", widehat(phi)[Gamma]), list(org=full.name1))}

  ylabel <- substitute(paste(org, " ", widehat(phi)), list(org=full.name2))
  if(gamma){ylabel <- substitute(paste(org, " ", widehat(phi)[Gamma]), list(org=full.name2))}

  plot(data[[org1.idx]], data[[org2.idx]],
       xlab = xlabel, ylab = ylabel)
    
  formula.reg <- ifelse(gamma,
                        paste(organism2, "_gamma_phi~", organism1, "_gamma_phi", sep=""),
                        paste(organism2, "_nogamma_phi~", organism1, "_nogamma_phi", sep=""))
  
  sma.model <- lmodel2(formula.reg, data = data, range.y=NULL, range.x=NULL, nperm=0)
  line.res <- as.numeric(sma.model$regression.results[3,2:3])
  abline(a = line.res[1], b = line.res[2])
  legend("bottomright", bty="n", legend = substitute(paste(rho, "=", r), list(r=round(sma.model$r, digits = 3))))
}


pdf("selac_phi_across_species_nogamma.pdf", width=10, height=10)
par(mfrow=c(5,5), mar=c(4,5,3,1)+0.1)

# change boolean flag to switch between gamma (TRUE) and no-gamma (FALSE) estimates of phi
gamma <- F

plotPhiCompSpecies(selac.exp.phi, "Scas", "Scer", "S. castellii", "S. cereviciae", gamma)
plotPhiCompSpecies(selac.exp.phi, "Skud", "Scer", "S. kudriavzevii", "S. cereviciae", gamma)
plotPhiCompSpecies(selac.exp.phi, "Smik", "Scer", "S. mikatae", "S. cereviciae", gamma)
plotPhiCompSpecies(selac.exp.phi, "Spar", "Scer", "S. paradoxus", "S. cereviciae", gamma)
plotPhiCompSpecies(selac.exp.phi, "Cgla", "Scer", "C. glabrata", "S. cereviciae", gamma)

plot(NULL, NULL)
plotPhiCompSpecies(selac.exp.phi, "Skud", "Scas", "S. kudriavzevii", "S. castellii", gamma)
plotPhiCompSpecies(selac.exp.phi, "Smik", "Scas", "S. mikatae", "S. castellii", gamma)
plotPhiCompSpecies(selac.exp.phi, "Spar", "Scas", "S. paradoxus", "S. castellii", gamma)
plotPhiCompSpecies(selac.exp.phi, "Cgla", "Scas", "C. glabrata", "S. castellii", gamma)

plot(NULL, NULL)
plot(NULL, NULL)
plotPhiCompSpecies(selac.exp.phi, "Smik", "Skud", "S. mikatae", "S. kudriavzevii", gamma)
plotPhiCompSpecies(selac.exp.phi, "Spar", "Skud", "S. paradoxus", "S. kudriavzevii", gamma)
plotPhiCompSpecies(selac.exp.phi, "Cgla", "Skud", "C. glabrata", "S. kudriavzevii", gamma)

plot(NULL, NULL)
plot(NULL, NULL)
plot(NULL, NULL)
plotPhiCompSpecies(selac.exp.phi, "Spar", "Smik", "S. paradoxus", "S. mikatae", gamma)
plotPhiCompSpecies(selac.exp.phi, "Cgla", "Smik", "C. glabrata", "S. mikatae", gamma)

plot(NULL, NULL)
plot(NULL, NULL)
plot(NULL, NULL)
plot(NULL, NULL)
plotPhiCompSpecies(selac.exp.phi, "Cgla", "Spar", "C. glabrata", "S. paradoxus", gamma)
dev.off()

######################
# COMAPRE EMPIRICAL EXPRESSION ACROSS SPECIES
######################


plotMeassureCompSpecies <- function(data, organism1, organism2, full.name1, full.name2, RNA.seq1, RNA.seq2){
  
  org1.idx <- ifelse(RNA.seq1, 
                     paste(organism1, "_RNA", sep=""),
                     paste(organism1, "_Microarray", sep=""))
  org2.idx <- ifelse(RNA.seq2, 
                     paste(organism2, "_RNA", sep=""),
                     paste(organism2, "_Microarray", sep=""))
  
  xlabel <- substitute(paste(log[10], " ", org[Microarray]), list(org=full.name1))
  if(RNA.seq1){xlabel <- substitute(paste(log[10], " ", org[RNA-Seq]), list(org=full.name1))}
  
  ylabel <- substitute(paste(log[10], " ", org[Microarray]), list(org=full.name2))
  if(RNA.seq2){ylabel <- substitute(paste(log[10], " ", org[RNA-Seq]), list(org=full.name2))}

  plot(log10(data[[org1.idx]]), log10(data[[org2.idx]]),
       xlab = xlabel, ylab = ylabel)
  
  formula.reg <- paste(organism2, ifelse(RNA.seq2, "_RNA", "_Microarray"), "~", 
                       organism1, ifelse(RNA.seq1, "_RNA", "_Microarray"), sep="")

  sma.model <- lmodel2(formula.reg, data = log10(data), range.y=NULL, range.x=NULL, nperm=0)
  line.res <- as.numeric(sma.model$regression.results[3,2:3])
  abline(a = line.res[1], b = line.res[2])
  legend("bottomright", bty="n", legend = substitute(paste(rho, "=", r), list(r=round(sma.model$r, digits = 3))))
}

pdf("selac_emp_comp.pdf", width=8, height=8)
par(mfrow=c(4,4), mar=c(4,5,3,1)+0.1)

plotMeassureCompSpecies(selac.exp.dat, "Scas", "Scer", "S. castellii", "S. cereviciae", F, T)
plotMeassureCompSpecies(selac.exp.dat, "Smik", "Scer", "S. mikatae", "S. cereviciae", T, T)
plotMeassureCompSpecies(selac.exp.dat, "Spar", "Scer", "S. paradoxus", "S. cereviciae", T, T)
plotMeassureCompSpecies(selac.exp.dat, "Cgla", "Scer", "C. glabrata", "S. cereviciae", F, T)

plot(NULL, NULL)
plotMeassureCompSpecies(selac.exp.dat, "Smik", "Scas", "S. mikatae", "S. castellii", T, F)
plotMeassureCompSpecies(selac.exp.dat, "Spar", "Scas", "S. paradoxus", "S. castellii", T, F)
plotMeassureCompSpecies(selac.exp.dat, "Cgla", "Scas", "C. glabrata", "S. castellii", F, F)

plot(NULL, NULL)
plot(NULL, NULL)
plotMeassureCompSpecies(selac.exp.dat, "Spar", "Smik", "S. paradoxus", "S. mikatae", T, T)
plotMeassureCompSpecies(selac.exp.dat, "Cgla", "Smik", "C. glabrata", "S. mikatae", F, T)

plot(NULL, NULL)
plot(NULL, NULL)
plot(NULL, NULL)
plotMeassureCompSpecies(selac.exp.dat, "Cgla", "Spar", "C. glabrata", "S. paradoxus", F, T)
dev.off()

######################
# COMAPRE EMPIRICAL EXPRESSION TO ESTIMATES
######################

plotEstVsEmp <- function(data.emp, data.est, full.name, desc.emp, gamma){
  ylabel <- substitute(paste(log[10], " ", org[d]), list(org=full.name, d=desc.emp))

  xlabel <- substitute(paste(log[10], " ", org, " ", phi), list(org=full.name))
  if(gamma){xlabel <- substitute(paste(log[10], " ", org, " ", phi[Gamma]), list(org=full.name))}
  
  plot(log10(data.est), log10(data.emp),
       xlab = xlabel, ylab = ylabel)
  
  data <- data.frame(emp=data.emp, est=data.est)

  sma.model <- lmodel2(emp~est, data = log10(data), range.y=NULL, range.x=NULL, nperm=0)
  line.res <- as.numeric(sma.model$regression.results[3,2:3])
  abline(a = line.res[1], b = line.res[2])
  legend("bottomright", bty="n", legend = substitute(paste(rho, "=", r), list(r=round(sma.model$r, digits = 3))))
  legend("topleft", bty="n", legend = full.name, text.font=3)
}

pdf("selac_phi_v_emp_nogamma.pdf", width=7, height=10)
par(mfrow=c(3,2), mar=c(4,5,3,1)+0.1)
plotEstVsEmp(selac.exp.dat$Scer_RNA, selac.exp.phi$Scer_nogamma_phi, "S. cereviciae", "RNA-Seq", T)
plotEstVsEmp(selac.exp.dat$Scas_Microarray, selac.exp.phi$Scas_nogamma_phi, "S. castellii", "Microarray", T)
plotEstVsEmp(selac.exp.dat$Smik_RNA, selac.exp.phi$Smik_nogamma_phi, "S. mikatae", "RNA-Seq", T)
plotEstVsEmp(selac.exp.dat$Spar_RNA, selac.exp.phi$Spar_nogamma_phi, "S. paradoxus", "RNA-Seq", T)
plotEstVsEmp(selac.exp.dat$Cgla_Microarray, selac.exp.phi$Cgla_nogamma_phi, "C. glabrata", "Microarray", T)
dev.off()

# change phi to psi in above function (plotEstVsEmp) for corect axis label
pdf("selac_psi_v_emp_gamma.pdf", width=7, height=10)
par(mfrow=c(3,2), mar=c(4,5,3,1)+0.1)
plotEstVsEmp(selac.exp.dat$Scer_RNA, selac.exp.dat$Psi_gamma, "S. cereviciae", "RNA-Seq", T)
plotEstVsEmp(selac.exp.dat$Scas_Microarray, selac.exp.dat$Psi_gamma, "S. castellii", "Microarray", T)
plotEstVsEmp(selac.exp.dat$Smik_RNA, selac.exp.dat$Psi_gamma, "S. mikatae", "RNA-Seq", T)
plotEstVsEmp(selac.exp.dat$Spar_RNA, selac.exp.dat$Psi_gamma, "S. paradoxus", "RNA-Seq", T)
plotEstVsEmp(selac.exp.dat$Cgla_Microarray, selac.exp.dat$Psi_gamma, "C. glabrata", "Microarray", T)
dev.off()
