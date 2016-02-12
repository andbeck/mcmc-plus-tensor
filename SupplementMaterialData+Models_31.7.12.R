# Requires the following libraries
library(mvtnorm)
library(MASS)
library(MCMCglmm)

##------------------------------------------------------------------------------------------------------
## Generate G matrices for multiple traits
##------------------------------------------------------------------------------------------------------
set.seed(123)
# All relationships highly positive between five traits measured 
# across two environments (10 x 10) (No Var Change No Cor Change)

G_1<-matrix(c(
  1.0,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,
  0.8,   1.0,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,
  0.8,   0.8,   1.0,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,
  0.8,   0.8,   0.8,   1.0,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,
  0.8,   0.8,   0.8,   0.8,   1.0,   0.8,   0.8,   0.8,   0.8,   0.8,
  0.8,   0.8,   0.8,   0.8,   0.8,   1.0,   0.8,   0.8,   0.8,   0.8,
  0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   1.0,   0.8,   0.8,   0.8,
  0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   1.0,   0.8,   0.8,
  0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   1.0,   0.8,
  0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   1.0),10,10) 

bv_1 <- lapply(1, mvrnorm, n=200, mu=rep(0,10), Sigma=G_1)
sire1<-lapply(bv_1, function(x) data.frame(sire=rep(1:nrow(x),each=4), rpt=rep(1:4,nrow(x)),
res1=rnorm(800,0,1.35^2),res2=rnorm(800,0,1.35^2),res3=rnorm(800,0,1.35^2),res4=rnorm(800,0,1.35^2),res5=rnorm(800,0,1.35^2),res6=rnorm(800,0,1.35^2),res7=rnorm(800,0,1.35^2),res8=rnorm(800,0,1.35^2),res9=rnorm(800,0,1.35^2),res10=rnorm(800,0,1.35^2), 
bv1=rep(x[,1],each=4),bv2=rep(x[,2],each=4),bv3=rep(x[,3],each=4),bv4=rep(x[,4],each=4),bv5=rep(x[,5],each=4),bv6=rep(x[,6],each=4),bv7=rep(x[,7],each=4),bv8=rep(x[,8],each=4),bv9=rep(x[,9],each=4),bv10=rep(x[,10],each=4)))
data1<-lapply(sire1, function(x) cbind(x, trait1=x$res1+x$bv1, trait2=x$res2+x$bv2, trait3=x$res3+x$bv3,trait4=x$res4+x$bv4,trait5=x$res5+x$bv5,trait6=x$res6+x$bv6,trait7=x$res7+x$bv7,trait8=x$res8+x$bv8,trait9=x$res9+x$bv9,trait10=x$res10+x$bv10))


# Relationships between traits constant within an environment but 
# are negative across environments (No Var Change - Cor Change)

G_2<-matrix(c(
  1.0,   0.8,   0.8,   0.8,   0.8,  -0.8,  -0.8,  -0.8,  -0.8,  -0.8,
  0.8,   1.0,   0.8,   0.8,   0.8,  -0.8,  -0.8,  -0.8,  -0.8,  -0.8,
  0.8,   0.8,   1.0,   0.8,   0.8,  -0.8,  -0.8,  -0.8,  -0.8,  -0.8,
  0.8,   0.8,   0.8,   1.0,   0.8,  -0.8,  -0.8,  -0.8,  -0.8,  -0.8,
  0.8,   0.8,   0.8,   0.8,   1.0,  -0.8,  -0.8,  -0.8,  -0.8,  -0.8,
 -0.8,  -0.8,  -0.8,  -0.8,  -0.8,   1.0,   0.8,   0.8,   0.8,   0.8,
 -0.8,  -0.8,  -0.8,  -0.8,  -0.8,   0.8,   1.0,   0.8,   0.8,   0.8,
 -0.8,  -0.8,  -0.8,  -0.8,  -0.8,   0.8,   0.8,   1.0,   0.8,   0.8,
 -0.8,  -0.8,  -0.8,  -0.8,  -0.8,   0.8,   0.8,   0.8,   1.0,   0.8,  
 -0.8,  -0.8,  -0.8,  -0.8,  -0.8,   0.8,   0.8,   0.8,   0.8,   1.0),10,10) 

bv_2 <- lapply(1, mvrnorm, n=200, mu=rep(0,10), Sigma=G_2)
sire2<-lapply(bv_2, function(x) data.frame(sire=rep(1:nrow(x),each=4), rpt=rep(1:4,nrow(x)),
res1=rnorm(800,0,1.35^2),res2=rnorm(800,0,1.35^2),res3=rnorm(800,0,1.35^2),res4=rnorm(800,0,1.35^2),res5=rnorm(800,0,1.35^2),res6=rnorm(800,0,1.35^2),res7=rnorm(800,0,1.35^2),res8=rnorm(800,0,1.35^2),res9=rnorm(800,0,1.35^2),res10=rnorm(800,0,1.35^2), 
bv1=rep(x[,1],each=4),bv2=rep(x[,2],each=4),bv3=rep(x[,3],each=4),bv4=rep(x[,4],each=4),bv5=rep(x[,5],each=4),bv6=rep(x[,6],each=4),bv7=rep(x[,7],each=4),bv8=rep(x[,8],each=4),bv9=rep(x[,9],each=4),bv10=rep(x[,10],each=4)))
data2<-lapply(sire2, function(x) cbind(x, trait1=x$res1+x$bv1, trait2=x$res2+x$bv2, trait3=x$res3+x$bv3,trait4=x$res4+x$bv4,trait5=x$res5+x$bv5,trait6=x$res6+x$bv6,trait7=x$res7+x$bv7,trait8=x$res8+x$bv8,trait9=x$res9+x$bv9,trait10=x$res10+x$bv10))

# Var Change, No Cor Change

G_3<-matrix(c(
	1.000, 0.800, 0.800, 0.800, 0.800, 0.800, 0.980, 1.131, 1.265, 1.386,
	0.800, 1.000, 0.800, 0.800, 0.800, 0.800, 0.980, 1.131, 1.265, 1.386,
	0.800, 0.800, 1.000, 0.800, 0.800, 0.800, 0.980, 1.131, 1.265, 1.386,
	0.800, 0.800, 0.800, 1.000, 0.800, 0.800, 0.980, 1.131, 1.265, 1.386,
	0.800, 0.800, 0.800, 0.800, 1.000, 0.800, 0.980, 1.131, 1.265, 1.386,
	0.800, 0.800, 0.800, 0.800, 0.800, 1.000, 0.980, 1.131, 1.265, 1.386,
	0.980, 0.980, 0.980, 0.980, 0.980, 0.980, 1.500, 1.386, 1.549, 1.697,
	1.131, 1.131, 1.131, 1.131, 1.131, 1.131, 1.386, 2.000, 1.789, 1.960,
	1.265, 1.265, 1.265, 1.265, 1.265, 1.265, 1.549, 1.789, 2.500, 2.191,
	1.386, 1.386, 1.386, 1.386, 1.386, 1.386, 1.697, 1.960, 2.191, 3.000),10,10)

bv_3 <- lapply(1, mvrnorm, n=200, mu=rep(0,10), Sigma=G_3)
sire3<-lapply(bv_3, function(x) data.frame(sire=rep(1:nrow(x),each=4), rpt=rep(1:4,nrow(x)),
res1=rnorm(800,0,1.35^2),res2=rnorm(800,0,1.35^2),res3=rnorm(800,0,1.35^2),res4=rnorm(800,0,1.35^2),res5=rnorm(800,0,1.35^2),res6=rnorm(800,0,1.35^2),res7=rnorm(800,0,1.35^2),res8=rnorm(800,0,1.35^2),res9=rnorm(800,0,1.35^2),res10=rnorm(800,0,1.35^2), 
bv1=rep(x[,1],each=4),bv2=rep(x[,2],each=4),bv3=rep(x[,3],each=4),bv4=rep(x[,4],each=4),bv5=rep(x[,5],each=4),bv6=rep(x[,6],each=4),bv7=rep(x[,7],each=4),bv8=rep(x[,8],each=4),bv9=rep(x[,9],each=4),bv10=rep(x[,10],each=4)))
data3<-lapply(sire3, function(x) cbind(x, trait1=x$res1+x$bv1, trait2=x$res2+x$bv2, trait3=x$res3+x$bv3,trait4=x$res4+x$bv4,trait5=x$res5+x$bv5,trait6=x$res6+x$bv6,trait7=x$res7+x$bv7,trait8=x$res8+x$bv8,trait9=x$res9+x$bv9,trait10=x$res10+x$bv10))

# One vector in environment 1, two vectors in environments 2 
# and relationships between traits reduce. Higher variation
# in environment 2 (Var Change and Cor Change)

G_4<-matrix(c(
  1.0,   0.8,   0.8,   0.8,   0.8,   0.8,   0.6,   0.4,   0.2,   0.0,
  0.8,   1.0,   0.8,   0.8,   0.8,   0.8,   0.6,   0.4,   0.2,   0.0,
  0.8,   0.8,   1.0,   0.8,   0.8,   0.8,   0.6,   0.4,   0.2,   0.0,
  0.8,   0.8,   0.8,   1.0,   0.8,   0.8,   0.6,   0.4,   0.2,   0.0,
  0.8,   0.8,   0.8,   0.8,   1.0,   0.8,   0.6,   0.4,   0.2,   0.0,
  0.8,   0.8,   0.8,   0.8,   0.8,   1.0,   0.6,   0.4,   0.2,   0.0,
  0.6,   0.6,   0.6,   0.6,   0.6,   0.6,   1.5,   0.8,   0.4,   0.2,
  0.4,   0.4,   0.4,   0.4,   0.4,   0.4,   0.8,   2.0,   1.0,   0.4,
  0.2,   0.2,   0.2,   0.2,   0.2,   0.2,   0.4,   1.0,   2.5,   1.2,
  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.2,   0.4,   1.2,   3.0),10,10) 

bv_4 <- lapply(1, mvrnorm, n=200, mu=rep(0,10), Sigma=G_4)
sire4<-lapply(bv_4, function(x) data.frame(sire=rep(1:nrow(x),each=4), rpt=rep(1:4,nrow(x)),
res1=rnorm(800,0,1.35^2),res2=rnorm(800,0,1.35^2),res3=rnorm(800,0,1.35^2),res4=rnorm(800,0,1.35^2),res5=rnorm(800,0,1.35^2),res6=rnorm(800,0,1.35^2),res7=rnorm(800,0,1.35^2),res8=rnorm(800,0,1.35^2),res9=rnorm(800,0,1.35^2),res10=rnorm(800,0,1.35^2), 
bv1=rep(x[,1],each=4),bv2=rep(x[,2],each=4),bv3=rep(x[,3],each=4),bv4=rep(x[,4],each=4),bv5=rep(x[,5],each=4),bv6=rep(x[,6],each=4),bv7=rep(x[,7],each=4),bv8=rep(x[,8],each=4),bv9=rep(x[,9],each=4),bv10=rep(x[,10],each=4)))
data4<-lapply(sire4, function(x) cbind(x, trait1=x$res1+x$bv1, trait2=x$res2+x$bv2, trait3=x$res3+x$bv3,trait4=x$res4+x$bv4,trait5=x$res5+x$bv5,trait6=x$res6+x$bv6,trait7=x$res7+x$bv7,trait8=x$res8+x$bv8,trait9=x$res9+x$bv9,trait10=x$res10+x$bv10))


# Two vectors in environment 1, two vectors in environments 2 
# and relationships between traits reduce. Higher variation
# in environment 2 

G_5<-matrix(c(
  1.0,   0.6,   0.4,   0.2,   0.0,  -0.8,  -0.6,  -0.4,  -0.2,   0.0,
  0.6,   1.0,   0.6,   0.4,   0.2,  -0.6,  -0.8,  -0.6,  -0.4,  -0.2,
  0.4,   0.6,   1.0,   0.6,   0.4,  -0.4,  -0.6,  -0.8,  -0.6,  -0.4,
  0.2,   0.4,   0.6,   1.0,   0.6,  -0.2,  -0.4,  -0.6,  -0.8,  -0.6,
  0.0,   0.2,   0.4,   0.6,   1.0,   0.0,  -0.2,  -0.4,  -0.6,  -0.8,
 -0.8,  -0.6,  -0.4,  -0.2,   0.0,   1.0,   0.6,   0.4,   0.2,   0.0,
 -0.6,  -0.8,  -0.6,  -0.4,  -0.2,   0.6,   1.5,   1.0,   0.4,   0.2,
 -0.4,  -0.6,  -0.8,  -0.6,  -0.4,   0.4,   1.0,   2.0,   1.2,   0.4,
 -0.2,  -0.4,  -0.6,  -0.8,  -0.6,   0.2,   0.4,   1.2,   2.5,   1.4,
  0.0,  -0.2,  -0.4,  -0.6,  -0.8,   0.0,   0.2,   0.4,   1.4,   3.0),10,10)

bv_5 <- lapply(1, mvrnorm, n=200, mu=rep(0,10), Sigma=G_5)
sire5<-lapply(bv_5, function(x) data.frame(sire=rep(1:nrow(x),each=4), rpt=rep(1:4,nrow(x)),
res1=rnorm(800,0,1.35^2),res2=rnorm(800,0,1.35^2),res3=rnorm(800,0,1.35^2),res4=rnorm(800,0,1.35^2),res5=rnorm(800,0,1.35^2),res6=rnorm(800,0,1.35^2),res7=rnorm(800,0,1.35^2),res8=rnorm(800,0,1.35^2),res9=rnorm(800,0,1.35^2),res10=rnorm(800,0,1.35^2), 
bv1=rep(x[,1],each=4),bv2=rep(x[,2],each=4),bv3=rep(x[,3],each=4),bv4=rep(x[,4],each=4),bv5=rep(x[,5],each=4),bv6=rep(x[,6],each=4),bv7=rep(x[,7],each=4),bv8=rep(x[,8],each=4),bv9=rep(x[,9],each=4),bv10=rep(x[,10],each=4)))
data5<-lapply(sire5, function(x) cbind(x, trait1=x$res1+x$bv1, trait2=x$res2+x$bv2, trait3=x$res3+x$bv3,trait4=x$res4+x$bv4,trait5=x$res5+x$bv5,trait6=x$res6+x$bv6,trait7=x$res7+x$bv7,trait8=x$res8+x$bv8,trait9=x$res9+x$bv9,trait10=x$res10+x$bv10))

##------------------------------------------------------------------------------------------------------
## Run Models for each Environment Using Multicore and get derived stats
##------------------------------------------------------------------------------------------------------

# set parameter expanded prior
prior1<-list(R=list(V=diag(5)*0.02, nu=6), 
			G=list(G1=list(V=diag(5),
						nu=5,
						alpha.mu=rep(0,5),
						alpha.V=diag(5)*1000)))

# load multicore package, rgl, derived stats function and psb plot function
library(multicore)
library(rgl)
source('~/Dropbox/###Subspace/Manuscript/AnalysisFunctions/derived.stats2.R', chdir = TRUE)
source('~/Dropbox/###Subspace/Manuscript/AnalysisFunctions/psb2.R', chdir = TRUE)

# all data made, and now collected into a list with indices.
alldat<-list(data1,data2,data3,data4,data5)

# create functions for E1 and E2
myMCMC.E1 <- function(i){
	use <- alldat[[i]][[1]]
	mod<-MCMCglmm(cbind(trait1, trait2, trait3, trait4, trait5)  ~ trait-1,
	random=~us(trait):sire, rcov=~us(trait):units,
	family=c("gaussian","gaussian","gaussian","gaussian", "gaussian"),
	prior=prior1, data=use, verbose=FALSE, burnin=50000, thin=100, nitt=150000)
	return(mod)
	}

myMCMC.E2 <- function(i){
	use <- alldat[[i]][[1]]
	mod<-MCMCglmm(cbind(trait6, trait7, trait8, trait9, trait10)  ~ trait-1,
	random=~us(trait):sire, rcov=~us(trait):units,
	family=c("gaussian","gaussian","gaussian","gaussian", "gaussian"),
	prior=prior1, data=use, verbose=FALSE, burnin=50000, thin=100, nitt=150000)
	return(mod)
	}

# Use mclapply to run on X cores - collect as list in Out.E1/2
Out.E1 <- mclapply(1:5, myMCMC.E1, mc.cores=5) # returns list
Out.E2 <- mclapply(1:5, myMCMC.E2, mc.cores=5) # returns list

# collect models into a list
mods<-list(Out.E1,Out.E2)

# function for derived stats on the mods list
myDerived<-function(i){
	d1<-mods[[1]][[i]]
	d2<-mods[[2]][[i]]
	return(derived.stats2(d1,d2, no.traits=5))
}

# collect set of derived stats
dstats<-mclapply(1:5, myDerived, mc.cores=5)

# Get the pictures
psb2(mods[[1]][[1]],mods[[2]][[1]], no.traits=5, method="G",shadeCA=FALSE)
psb2(mods[[1]][[2]],mods[[2]][[2]], no.traits=5, method="G",shadeCA=FALSE)
psb2(mods[[1]][[3]],mods[[2]][[3]], no.traits=5, method="G",shadeCB=FALSE)
psb2(mods[[1]][[4]],mods[[2]][[4]], no.traits=5, method="G",shadeCB=FALSE)
psb2(mods[[1]][[5]],mods[[2]][[5]], no.traits=5, method="G",shadeCB=FALSE)

##------------------------------------------------------------------------------------------------------
## The long way
##------------------------------------------------------------------------------------------------------

## models
model1.1<-MCMCglmm(cbind(trait1, trait2, trait3, trait4, trait5)  ~ trait-1, random=~us(trait):sire, rcov=~us(trait):units, family=c("gaussian","gaussian","gaussian","gaussian", "gaussian"), prior=prior1, data=data1[[1]], verbose=FALSE, burnin=50000, thin=100, nitt=150000)
model1.2<-MCMCglmm(cbind(trait6, trait7, trait8, trait9, trait10) ~ trait-1, random=~us(trait):sire, rcov=~us(trait):units, family=c("gaussian","gaussian","gaussian","gaussian", "gaussian"), prior=prior1, data=data1[[1]], verbose=FALSE, burnin=50000, thin=100, nitt=150000)

model2.1<-MCMCglmm(cbind(trait1, trait2, trait3, trait4, trait5)  ~ trait-1, random=~us(trait):sire, rcov=~us(trait):units, family=c("gaussian","gaussian","gaussian","gaussian", "gaussian"), prior=prior1, data=data2[[1]], verbose=FALSE, burnin=50000, thin=100, nitt=150000)
model2.2<-MCMCglmm(cbind(trait6, trait7, trait8, trait9, trait10) ~ trait-1, random=~us(trait):sire, rcov=~us(trait):units, family=c("gaussian","gaussian","gaussian","gaussian", "gaussian"), prior=prior1, data=data2[[1]], verbose=FALSE, burnin=50000, thin=100, nitt=150000)

model3.1<-MCMCglmm(cbind(trait1, trait2, trait3, trait4, trait5)  ~ trait-1, random=~us(trait):sire, rcov=~us(trait):units, family=c("gaussian","gaussian","gaussian","gaussian", "gaussian"), prior=prior1, data=data3[[1]], verbose=FALSE, burnin=50000, thin=100, nitt=150000)
model3.2<-MCMCglmm(cbind(trait6, trait7, trait8, trait9, trait10) ~ trait-1, random=~us(trait):sire, rcov=~us(trait):units, family=c("gaussian","gaussian","gaussian","gaussian", "gaussian"), prior=prior1, data=data3[[1]], verbose=FALSE, burnin=50000, thin=100, nitt=150000)

model4.1<-MCMCglmm(cbind(trait1, trait2, trait3, trait4, trait5)  ~ trait-1, random=~us(trait):sire, rcov=~us(trait):units, family=c("gaussian","gaussian","gaussian","gaussian", "gaussian"), prior=prior1, data=data4[[1]], verbose=FALSE, burnin=50000, thin=100, nitt=150000)
model4.2<-MCMCglmm(cbind(trait6, trait7, trait8, trait9, trait10) ~ trait-1, random=~us(trait):sire, rcov=~us(trait):units, family=c("gaussian","gaussian","gaussian","gaussian", "gaussian"), prior=prior1, data=data4[[1]], verbose=FALSE, burnin=50000, thin=100, nitt=150000)

model5.1<-MCMCglmm(cbind(trait1, trait2, trait3, trait4, trait5)  ~ trait-1, random=~us(trait):sire, rcov=~us(trait):units, family=c("gaussian","gaussian","gaussian","gaussian", "gaussian"), prior=prior1, data=data4[[1]], verbose=FALSE, burnin=50000, thin=100, nitt=150000)
model5.2<-MCMCglmm(cbind(trait6, trait7, trait8, trait9, trait10) ~ trait-1, random=~us(trait):sire, rcov=~us(trait):units, family=c("gaussian","gaussian","gaussian","gaussian", "gaussian"), prior=prior1, data=data4[[1]], verbose=FALSE, burnin=50000, thin=100, nitt=150000)

drs1<-derived.stats2(model1.1,model1.2, no.traits=5)
drs2<-derived.stats2(model2.1,model2.2, no.traits=5)
drs3<-derived.stats2(model3.1,model3.2, no.traits=5)
drs4<-derived.stats2(model4.1,model4.2, no.traits=5)