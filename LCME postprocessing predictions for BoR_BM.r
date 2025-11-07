### LCME post-processing predictions ###
### 1. Get LCME covariates
### 2. Resample LCME posteriors
### 3. Predict population using LCME posteriors and covariates

#setwd('D:/FWS/JAGS/CAA/DSLCM/Models/DSLCM_Entrainment/Output/')
setwd("~/GitHub/DeltaSmelt_LCM")
load('Dat_LCM3_run_v4_2020_11_23.RData') # best model
library(R2jags) # For jags.parallel(). Will also load required packages rjags, coda, and link to JAGS software
library(abind)
out2<-as.mcmc(fitted.model)

n.samp<-1e5 # number of simulations and posterior samples

wtr.yr<-c(1,1,1,1,1,1,4,4,2,3,3,1,4,5,4,3,1,3,4,5,5) # Sacto WY type wet = 1, critical = 5

## Part 1: covariates

# Action effects
#CalSim <- read.csv('D:/FWS/JAGS/CAA/DSLCM/Models/DSLCM_Entrainment/BoR Reiniation 2023/FlowZoopData_2022ROC_EffectsAnalysis_CohortYear.csv',header=T)
CalSim <- read.csv(file.path("Salinity_Zooplankton_analysis","FlowZoopData_2022ROC_EffectsAnalysis_CohortYear_2024-09-25.csv"),header=T)

CalSim[,5] <- 1.9834710990151385*CalSim[,5] # convert cfs to sum of acre-feet for 3 months
CalSim.act <- c('Alt1','Alt2v1wTUCP','Alt2v1woTUCP','Alt2v2noTUCP','Alt2v3noTUCP','Alt3','Alt4','Alt5','EXP1','EXP3','NAA')

# Covariates
x1<-cbind(rec2.cov.mat.stand,PL2.m.cov.mat.stand,Juv.m.cov.mat.stand,PL1.f.cov.mat.stand[,(1:2)],PL2.f.cov.mat.stand[,(1:2)])
x1.unstand<-cbind(rec2.cov.mat,PL2.m.cov.mat,Juv.m.cov.mat,PL1.f.cov.mat[,(1:2)],PL2.f.cov.mat[,(1:2)])
x2<-cbind(SA1.m.cov.mat.stand,SA2.m.cov.mat.stand,A1.m.cov.mat.stand)
for(i in 1:nrow(SA1.m.cov.mat)) {
	SA1.m.cov.mat[i,1] <- max(1,SA1.m.cov.mat[i,1]) # Change 0s to 1
	}
x2.unstand<-cbind(SA1.m.cov.mat,SA2.m.cov.mat,A1.m.cov.mat)
x3<-cbind(SA1.f.cov.mat.stand[,(1:2)],SA2.f.cov.mat.stand[,(1:2)],A1.f.cov.mat.stand[,(1:2)])
x3.unstand<-cbind(SA1.f.cov.mat[,(1:2)],SA2.f.cov.mat[,(1:2)],A1.f.cov.mat[,(1:2)])

XRec2mn <- mean(rec2.cov.mat)
XMPL2mn <- mean(PL2.m.cov.mat)
XMJmn <- mean(Juv.m.cov.mat)
XMSA1mn <- mean(SA1.m.cov.mat)
XMSA2mn <- mean(c(SA2.m.cov.mat,A1.m.cov.mat))
XMA1mn <- mean(c(SA2.m.cov.mat,A1.m.cov.mat))
XRec2sd <- sd(rec2.cov.mat)
XMPL2sd <- sd(PL2.m.cov.mat)
XMJsd <- sd(Juv.m.cov.mat)
XMSA1sd <- sd(SA1.m.cov.mat)
XMSA2sd <- sd(c(SA2.m.cov.mat,A1.m.cov.mat))
XMA1sd <- sd(c(SA2.m.cov.mat,A1.m.cov.mat))

XFPLmn <- c(mean(c(PL1.f.cov.mat[,1],PL2.f.cov.mat[,1])),mean(c(PL1.f.cov.mat[,2],PL2.f.cov.mat[,2])))
XFSAmn <- c(mean(c(SA1.f.cov.mat[,1],SA2.f.cov.mat[,1],A1.f.cov.mat[,1])),mean(c(SA1.f.cov.mat[,2],SA2.f.cov.mat[,2],A1.f.cov.mat[,2])))
XFPLsd <- c(sd(c(PL1.f.cov.mat[,1],PL2.f.cov.mat[,1])),sd(c(PL1.f.cov.mat[,2],PL1.f.cov.mat[,2])))
XFSAsd <- c(sd(c(SA1.f.cov.mat[,1],SA2.f.cov.mat[,1],A1.f.cov.mat[,1])),sd(c(SA1.f.cov.mat[,2],SA2.f.cov.mat[,2],A1.f.cov.mat[,2])))

x<-cbind(x1,x2,x3)
x.unstand<-cbind(x1.unstand,x2.unstand,x3.unstand)

new.SDT<-array(NA,dim=c(Y,5))
new.OMR<-array(NA,dim=c(Y,5))

### Part 2: sample LCME posteriors
out2<-as.mcmc(fitted.model)
set.seed(123)
rand.index<-sample((1:nrow(as.matrix(out2))),n.samp,replace=TRUE)
post.samp<-as.matrix(out2)[rand.index,]
Ninit<-post.samp[,'Ninit']
sig.repro<-post.samp[,'sig.repro[1]']
sig.M1<-post.samp[,'sig.M[1]']
sig.M4<-post.samp[,'sig.M[4]']
sig.F<-post.samp[,'sig.F[1]']
alphaR11<-post.samp[,'alphaR1[1]']
alphaR21<-post.samp[,'alphaR2[1]']
alphaR22<-post.samp[,'alphaR2[2]']
betaPL11<-post.samp[,'betaPL1[1]']
betaPL21<-post.samp[,'betaPL2[1]']
betaPL22<-post.samp[,'betaPL2[2]']
betaJ1<-post.samp[,'betaJ[1]']
betaJ2<-post.samp[,'betaJ[2]']
betaSA11<-post.samp[,'betaSA1[1]']
betaSA12<-post.samp[,'betaSA1[2]']
betaSA21<-post.samp[,'betaSA2[1]']
betaSA22<-post.samp[,'betaSA2[2]']
betaA11<-post.samp[,'betaA1[1]']
betaA12<-post.samp[,'betaA1[2]']
gammaPL11<-post.samp[,'gammaPL1[1]']
gammaPL12<-post.samp[,'gammaPL1[2]']
gammaPL13<-post.samp[,'gammaPL1[3]']
gammaPL14<-post.samp[,'gammaPL1[4]']
gammaPL21<-post.samp[,'gammaPL2[1]']
gammaPL22<-post.samp[,'gammaPL2[2]']
gammaPL23<-post.samp[,'gammaPL2[3]']
gammaPL24<-post.samp[,'gammaPL2[4]']
gammaSA11<-post.samp[,'gammaSA1[1]']
gammaSA12<-post.samp[,'gammaSA1[2]']
gammaSA13<-post.samp[,'gammaSA1[3]']
gammaSA21<-post.samp[,'gammaSA2[1]']
gammaSA22<-post.samp[,'gammaSA2[2]']
gammaSA23<-post.samp[,'gammaSA2[3]']
gammaA11<-post.samp[,'gammaA1[1]']
gammaA12<-post.samp[,'gammaA1[2]']
gammaA13<-post.samp[,'gammaA1[3]']
B1<-post.samp[,'B[1]']
B2<-post.samp[,'B[2]']

all.lambda<-all.lamAB<-all.lam.mn <- NULL

for (k in 1:length(CalSim.act)) { # loop through each action
XMPL2.act <- CalSim[which(CalSim[,2]==paste0(CalSim.act[k])),5]
XF.act <- cbind(CalSim[which(CalSim[,2]==paste0(CalSim.act[k])),6],
 CalSim[which(CalSim[,2]==paste0(CalSim.act[k])),7],
 CalSim[which(CalSim[,2]==paste0(CalSim.act[k])),8],
 CalSim[which(CalSim[,2]==paste0(CalSim.act[k])),9],
 CalSim[which(CalSim[,2]==paste0(CalSim.act[k])),10])
XMSA2.act <- CalSim[which(CalSim[,2]==paste0(CalSim.act[k])),3]
XMA1.act <- CalSim[which(CalSim[,2]==paste0(CalSim.act[k])),4]

new.XRec2<-x[,1]
#for (y in 1:Y) {
# new.XRec2[y]<-if(is.na(XRec2.act[y])) {x[y,1]} else {(XRec2.act[y]-XRec2mn)/XRec2sd}
# }
new.XMPL2<-x[,2]
for (y in 1:Y) {
 new.XMPL2[y]<-if(is.na(XMPL2.act[y])) {x[y,2]} else {(XMPL2.act[y]-XMPL2mn)/XMPL2sd}
 }
new.XMJ<-x[,3]
#for (y in 1:Y) {
# new.XMJ[y]<-if(is.na(XMJ.act[y])) {x[y,3]} else {(XMJ.act[y]-XMJmn)/XMJsd}
# }
new.SDT[,1]<-x[,4]
new.SDT[,2]<-x[,6]
new.SDT[,3]<-x[,11]
new.SDT[,4]<-x[,13]
new.SDT[,5]<-x[,15]
new.OMR[,1]<-x[,5]
new.OMR[,2]<-x[,7]
new.OMR[,3]<-x[,12]
new.OMR[,4]<-x[,14]
new.OMR[,5]<-x[,16]
for (y in 1:Y) {
 new.OMR[y,1]<-if(is.na(XF.act[y,1])) {x[y,5]} else {(XF.act[y,1]-XFPLmn[2])/XFPLsd[2]}
 new.OMR[y,2]<-if(is.na(XF.act[y,2])) {x[y,7]} else {(XF.act[y,2]-XFPLmn[2])/XFPLsd[2]}
 new.OMR[y,3]<-if(is.na(XF.act[y,3])) {x[y,12]} else {(XF.act[y,3]-XFSAmn[2])/XFSAsd[2]}
 new.OMR[y,4]<-if(is.na(XF.act[y,4])) {x[y,14]} else {(XF.act[y,4]-XFSAmn[2])/XFSAsd[2]}
 new.OMR[y,5]<-if(is.na(XF.act[y,5])) {x[y,16]} else {(XF.act[y,5]-XFSAmn[2])/XFSAsd[2]}
 }
new.XMSA1<-x[,8]
new.XMSA2<-x[,9]
for (y in 1:Y) {
 new.XMSA2[y]<-if(is.na(XMSA2.act[y])) {x[y,9]} else {(XMSA2.act[y]-XMSA2mn)/XMSA2sd}
 }
new.XMA1<-x[,10]
for (y in 1:Y) {
 new.XMA1[y]<-if(is.na(XMA1.act[y])) {x[y,10]} else {(XMA1.act[y]-XMA1mn)/XMA1sd}
 }

### Part 3: LCME predictions
new.lambda <- matrix(NA,Y,n.samp)

# ############# State-process model ##############
# Spawning and survival of initial A1 cohort
Ninit.SSB1 <- Ninit*Ninit.SW
Ninit.M <- rlnorm(n.samp,betaA11,sig.M4)
Ninit.F <- rlnorm(n.samp,gammaA11,sig.M4)

Ninit.phi <- exp(-(Ninit.M+Ninit.F))
Ninit.SSB2 <- Ninit*Ninit.phi*Apr.SW[1]

# Spawner-recruit model
# nPL1 and nPL2 given nA1 and nA2
for (y in 1:1) { # 1 variable models
  
 mu.rho1 <- alphaR11 #+inprod(rec.cov1[y,],alphaR1[2:(1+num_unique_cov_rec1)])
 mu.rho2 <- alphaR21+new.XRec2[y]*alphaR22

 rho1 <- rlnorm(n.samp,mu.rho1,sig.repro)
 rho2 <- rlnorm(n.samp,mu.rho2,sig.repro)

# survival (phi) model

# Entrainmnet mortality
 mu.F1 <- gammaPL11+new.SDT[y,1]*gammaPL12+new.OMR[y,1]*gammaPL13+new.SDT[y,1]*new.OMR[y,1]*gammaPL14
 mu.F2 <- gammaPL21+new.SDT[y,2]*gammaPL22+new.OMR[y,2]*gammaPL23+new.SDT[y,2]*new.OMR[y,2]*gammaPL24
 mu.F4 <- gammaSA11+new.SDT[y,3]*gammaSA12+new.OMR[y,3]*gammaSA13
 mu.F5 <- gammaSA21+new.SDT[y,4]*gammaSA12+new.OMR[y,4]*gammaSA13
 mu.F6 <- gammaA11+new.SDT[y,5]*gammaSA12+new.OMR[y,5]*gammaSA13
 F1 <- rlnorm(n.samp,mu.F1,sig.F)
 F2 <- rlnorm(n.samp,mu.F2,sig.F)
 F3 <- 0
 F4 <- rlnorm(n.samp,mu.F4,sig.F)
 F5 <- rlnorm(n.samp,mu.F5,sig.F)
 F6 <- rlnorm(n.samp,mu.F6,sig.F)
 
# Natural mortality
 mu.M1 <- betaPL11
 mu.M2 <- betaPL21+new.XMPL2[y]*betaPL22
 mu.M3 <- betaJ1+new.XMJ[y]*betaJ2
 mu.M4 <- betaSA11+new.XMSA1[y]*betaSA12
 mu.M5 <- betaSA21+new.XMSA2[y]*betaSA22
 mu.M6 <- betaA11+new.XMA1[y]*betaA12
  
 M1 <- rlnorm(n.samp,mu.M1,sig.M1)
 M2 <- rlnorm(n.samp,mu.M2,sig.M1)
 M3 <- rlnorm(n.samp,mu.M3,sig.M1)
 M4 <- rlnorm(n.samp,mu.M4,sig.M4)
 M5 <- rlnorm(n.samp,mu.M5,sig.M4)
 M6 <- rlnorm(n.samp,mu.M6,sig.M4)

# Survival
 Z1 <- F1+M1
 Z2 <- F2+M2
 Z3 <- F3+M3
 Z4 <- F4+M4
 Z5 <- F5+M5
 Z6 <- F6+M6
 
 phi1 <- exp(-Z1)
 phi2 <- exp(-Z2)
 phi3 <- exp(-Z3)
 phi4 <- exp(-Z4)
 phi5 <- exp(-Z5)
 phi6 <- exp(-Z6)

# Abundance
 nPL1 <- Ninit.SSB1*rho1
 nPL2 <- nPL1*phi1+Ninit.SSB2*rho2
 nJ <- nPL2*phi2
 nSA1 <- nJ*phi3
 nSA2 <- nSA1*phi4
 nA1 <- nSA2*phi5
 nA2 <- nA1*phi6
 
 SSB1 <- nA1*Mar.SW[y] # length-wt conversion
 SSB2 <- nA2*Apr.SW[y]
 
 new.lambda[1,] <- nA1/Ninit
 priornA1 <- nA1
 }

# Spawner-recruit model
for (y in 2:Y) { 
 mu.rho1 <- alphaR11 #+inprod(rec.cov1[y,],alphaR1[2:(1+num_unique_cov_rec1)])
 mu.rho2 <- alphaR21+new.XRec2[y]*alphaR22

 rho1 <- rlnorm(n.samp,mu.rho1,sig.repro)
 rho2 <- rlnorm(n.samp,mu.rho2,sig.repro)

# survival (phi) model

# Entrainmnet mortality
 mu.F1 <- gammaPL11+new.SDT[y,1]*gammaPL12+new.OMR[y,1]*gammaPL13+new.SDT[y,1]*new.OMR[y,1]*gammaPL14
 mu.F2 <- gammaPL21+new.SDT[y,2]*gammaPL22+new.OMR[y,2]*gammaPL23+new.SDT[y,2]*new.OMR[y,2]*gammaPL24
 mu.F4 <- gammaSA11+new.SDT[y,3]*gammaSA12+new.OMR[y,3]*gammaSA13
 mu.F5 <- gammaSA21+new.SDT[y,4]*gammaSA12+new.OMR[y,4]*gammaSA13
 mu.F6 <- gammaA11+new.SDT[y,5]*gammaSA12+new.OMR[y,5]*gammaSA13
 F1 <- rlnorm(n.samp,mu.F1,sig.F)
 F2 <- rlnorm(n.samp,mu.F2,sig.F)
 F3 <- 0
 F4 <- rlnorm(n.samp,mu.F4,sig.F)
 F5 <- rlnorm(n.samp,mu.F5,sig.F)
 F6 <- rlnorm(n.samp,mu.F6,sig.F)
 
# Natural mortality
 mu.M1 <- betaPL11
 mu.M2 <- betaPL21+new.XMPL2[y]*betaPL22
 mu.M3 <- betaJ1+new.XMJ[y]*betaJ2
 mu.M4 <- betaSA11+new.XMSA1[y]*betaSA12
 mu.M5 <- betaSA21+new.XMSA2[y]*betaSA22
 mu.M6 <- betaA11+new.XMA1[y]*betaA12
  
 M1 <- rlnorm(n.samp,mu.M1,sig.M1)
 M2 <- rlnorm(n.samp,mu.M2,sig.M1)
 M3 <- rlnorm(n.samp,mu.M3,sig.M1)
 M4 <- rlnorm(n.samp,mu.M4,sig.M4)
 M5 <- rlnorm(n.samp,mu.M5,sig.M4)
 M6 <- rlnorm(n.samp,mu.M6,sig.M4)

# Survival
 Z1 <- F1+M1
 Z2 <- F2+M2
 Z3 <- F3+M3
 Z4 <- F4+M4
 Z5 <- F5+M5
 Z6 <- F6+M6
 
 phi1 <- exp(-Z1)
 phi2 <- exp(-Z2)
 phi3 <- exp(-Z3)
 phi4 <- exp(-Z4)
 phi5 <- exp(-Z5)
 phi6 <- exp(-Z6)
 
# Abundance
 nPL1 <- SSB1*rho1
 nPL2 <- nPL1*phi1+SSB2*rho2
 nJ <- nPL2*phi2
 nSA1 <- nJ*phi3
 nSA2 <- nSA1*phi4
 nA1 <- nSA2*phi5
 nA2 <- nA1*phi6
 
 SSB1 <- nA1*Mar.SW[y] # length-wt conversion
 SSB2 <- nA2*Apr.SW[y]
 
 new.lambda[y,] <- nA1/priornA1
 priornA1 <- nA1
 }

new.lambda<-t(new.lambda)
yr.seq<-seq(cohort.analysis.range[1],cohort.analysis.range[2],by=1)
colnames(new.lambda)<-yr.seq

boxplot(new.lambda,outline=F)

lamAB<-matrix(NA,Y,11)
lam.mn<-vector()
for (t in 1:Y) { 
 lamAB[t,1] <- mean(new.lambda[,t],na.rm=T)
 lamAB[t,2] <- min(new.lambda[,t],na.rm=T)
 lamAB[t,3] <- max(new.lambda[,t],na.rm=T)
 lamAB[t,4] <- quantile(new.lambda[,t],0.025,na.rm=T)
 lamAB[t,5] <- quantile(new.lambda[,t],0.05,na.rm=T)
 lamAB[t,6] <- quantile(new.lambda[,t],0.25,na.rm=T)
 lamAB[t,7] <- quantile(new.lambda[,t],0.5,na.rm=T)
 lamAB[t,8] <- quantile(new.lambda[,t],0.75,na.rm=T)
 lamAB[t,9] <- quantile(new.lambda[,t],0.90,na.rm=T)
 lamAB[t,10] <- quantile(new.lambda[,t],0.975,na.rm=T)
 lamAB[t,11] <- (sd(new.lambda[,t],na.rm=T)/mean(new.lambda[,t],na.rm=T))
 }
colnames(lamAB) <- c("mean","min","max","2.5%","5%","25%","50%","75%","90%","97.5%","CV")

lam.mn[1] <- exp(mean(log(lamAB[,7]))) # print geometric mean pop. growth rate 1995-2015
lam.mn[2] <- exp(mean(log(lamAB[12:Y,7]))) # print geometric mean pop. growth rate 2007-2015
lam.mn[3] <- exp(mean(log(lamAB[11:Y,7]))) # print geometric mean pop. growth rate 2006-2015
lam.mn[4] <- exp(mean(log(lamAB[1:10,7]))) # print geometric mean pop. growth rate 1995-2005
lam.mn[5] <- exp(mean(log(lamAB[which(wtr.yr[1:Y]<=2),7]),na.rm=T)) # print geometric mean pop. growth rate AN and wet yrs
lam.mn[6] <- exp(mean(log(lamAB[which(wtr.yr[1:Y]>=4),7]),na.rm=T)) # print geometric mean pop. growth rate dry and critical yrs

all.lambda <- abind(all.lambda,new.lambda,along=3)
all.lamAB <- cbind(all.lamAB,lamAB[,7])
all.lam.mn <- cbind(all.lam.mn,lam.mn)
}

colnames(all.lamAB)<-colnames(all.lam.mn)<-CalSim.act
#write.table(all.lamAB,file='C:/Users/wsmith/Desktop/lamAB.ROC.csv')
#write.table(all.lam.mn,file='C:/Users/wsmith/Desktop/lam.mn.ROC.csv')
write.csv(all.lamAB,file=file.path("DeltaSmeltLCME_output","lamAB.ROC_2025-06-13.csv"))
write.csv(all.lam.mn,file=file.path("DeltaSmeltLCME_output","lam.mn.ROC_2025-06-13.csv"))

#write.table(all.lam.mn,file='lam.mn.ROC.csv')
