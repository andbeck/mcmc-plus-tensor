# load the libraries required
library(mvtnorm)
library(MASS)
library(MCMCglmm)
library(corpcor)

## Example 1: two traits with no covariance and no GEI in either trait across E
G_A<-matrix(c(
   1.0,  0.0,  0.0,  0.0,   
   0.0,  0.0,  0.0,  0.0, 
   0.0,  0.0,  1.0,  0.0,
   0.0,  0.0,  0.0,  0.0),4,4)

bv_a <- lapply(1, mvrnorm, n=200, mu=rep(0,4), Sigma=G_A)


# simulate a dataframe from the G matrix 
sire1<-lapply(bv_a, function(x) data.frame(ID=rep(1:nrow(x),each=20), time=rep(-2:2,nrow(x)), intercept1=rep(x[,1],each=20),slope1=rep(x[,2],each=20),
res1=rnorm(4000,0,1.35^2), intercept2=rep(x[,3],each=20),slope2=rep(x[,4],each=20),res2=rnorm(4000,0,1.35^2)))

data1<-lapply(sire1, function(x) cbind(ID=x$ID, time=x$time, phen1=x$intercept1+(x$slope1*x$time)+x$res1, phen2=x$intercept2+(x$slope2*x$time)+x$res2))
data1<-data.frame(data1[[1]])
data1$ID<-as.factor(data1$ID)
data1$group<-as.factor(data1$time)
data1<-data1[order(data1$ID,data1$time),]
head(data1)

d1<-data1[data1$time==-2,]
d2<-data1[data1$time==-1,]
d3<-data1[data1$time==-0,]
d4<-data1[data1$time== 1,]
d5<-data1[data1$time== 2,]

prior1<-list(R=list(V=diag(2)*0.02, nu=3), 
			G=list(G1=list(V=diag(2),
						nu=2,
						alpha.mu=rep(0,2),
						alpha.V=diag(2)*1000)))

m1a<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d1,
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m2a<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d2, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m3a<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d3, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m4a<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d4, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m5a<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d5, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)



## Example 2: two traits no genetic covariance between them, GEI for one trait but not the other
# simulate breeding values for 200 genotypes, for two traits. Each trait has intercept variance of 1; one trait has GEI with variance of 0.5; the other has no GEI variance; traits are uncorrelated
G_A<-matrix(c(
   1.0,  0.0,  0.0,  0.0,   
   0.0,  0.5,  0.0,  0.0, 
   0.0,  0.0,  1.0,  0.0,
   0.0,  0.0,  0.0,  0.0),4,4)

bv_a <- lapply(1, mvrnorm, n=200, mu=rep(0,4), Sigma=G_A)


# simulate a dataframe from the G matrix 
sire1<-lapply(bv_a, function(x) data.frame(ID=rep(1:nrow(x),each=20), time=rep(-2:2,nrow(x)), intercept1=rep(x[,1],each=20),slope1=rep(x[,2],each=20),
res1=rnorm(4000,0,1.35^2), intercept2=rep(x[,3],each=20),slope2=rep(x[,4],each=20),res2=rnorm(4000,0,1.35^2)))

data1<-lapply(sire1, function(x) cbind(ID=x$ID, time=x$time, phen1=x$intercept1+(x$slope1*x$time)+x$res1, phen2=x$intercept2+(x$slope2*x$time)+x$res2))
data1<-data.frame(data1[[1]])
data1$ID<-as.factor(data1$ID)
data1$group<-as.factor(data1$time)
data1<-data1[order(data1$ID,data1$time),]
head(data1)

d1<-data1[data1$time==-2,]
d2<-data1[data1$time==-1,]
d3<-data1[data1$time==-0,]
d4<-data1[data1$time== 1,]
d5<-data1[data1$time== 2,]

prior1<-list(R=list(V=diag(2)*0.02, nu=3), 
			G=list(G1=list(V=diag(2),
						nu=2,
						alpha.mu=rep(0,2),
						alpha.V=diag(2)*1000)))

m1b<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d1,
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m2b<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d2, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m3b<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d3, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m4b<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d4, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m5b<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d5, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

## Example 3. two traits, with zero genetic covariance, both with GEI
# simulate breeding values for 200 genotypes, for two traits. Each trait has intercept variance of 1; one trait has GEI with variance of 0.5; the other has no GEI variance; traits are uncorrelated
G_A<-matrix(c(
   1.0,  0.0,  0.0,  0.0,   
   0.0,  0.5,  0.0,  0.0, 
   0.0,  0.0,  1.0,  0.0,
   0.0,  0.0,  0.0,  0.5),4,4)

bv_a <- lapply(1, mvrnorm, n=200, mu=rep(0,4), Sigma=G_A)


# simulate a dataframe from the G matrix 
sire1<-lapply(bv_a, function(x) data.frame(ID=rep(1:nrow(x),each=20), time=rep(-2:2,nrow(x)), intercept1=rep(x[,1],each=20),slope1=rep(x[,2],each=20),
res1=rnorm(4000,0,1.35^2), intercept2=rep(x[,3],each=20),slope2=rep(x[,4],each=20),res2=rnorm(4000,0,1.35^2)))

data1<-lapply(sire1, function(x) cbind(ID=x$ID, time=x$time, phen1=x$intercept1+(x$slope1*x$time)+x$res1, phen2=x$intercept2+(x$slope2*x$time)+x$res2))
data1<-data.frame(data1[[1]])
data1$ID<-as.factor(data1$ID)
data1$group<-as.factor(data1$time)
data1<-data1[order(data1$ID,data1$time),]
head(data1)

d1<-data1[data1$time==-2,]
d2<-data1[data1$time==-1,]
d3<-data1[data1$time==-0,]
d4<-data1[data1$time== 1,]
d5<-data1[data1$time== 2,]

prior1<-list(R=list(V=diag(2)*0.02, nu=3), 
			G=list(G1=list(V=diag(2),
						nu=2,
						alpha.mu=rep(0,2),
						alpha.V=diag(2)*1000)))

m1c<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d1,
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m2c<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d2, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m3c<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d3, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m4c<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d4, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m5c<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d5, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

## Example 4. two traits, with zero genetic covariance, both with GEI
# simulate breeding values for 200 genotypes, for two traits. Each trait has intercept variance of 1; one trait has GEI with variance of 0.5; the other has no GEI variance; traits are uncorrelated
G_A<-matrix(c(
   1.0,  0.0,  0.9,  0.64,   
   0.0,  0.5,  0.64,  0.45, 
   0.9,  0.64,  1.0,  0.0,
   0.64,  0.45,  0.0,  0.5),4,4)

bv_a <- lapply(1, mvrnorm, n=200, mu=rep(0,4), Sigma=G_A)


# simulate a dataframe from the G matrix 
sire1<-lapply(bv_a, function(x) data.frame(ID=rep(1:nrow(x),each=20), time=rep(-2:2,nrow(x)), intercept1=rep(x[,1],each=20),slope1=rep(x[,2],each=20),
res1=rnorm(4000,0,1.35^2), intercept2=rep(x[,3],each=20),slope2=rep(x[,4],each=20),res2=rnorm(4000,0,1.35^2)))

data1<-lapply(sire1, function(x) cbind(ID=x$ID, time=x$time, phen1=x$intercept1+(x$slope1*x$time)+x$res1, phen2=x$intercept2+(x$slope2*x$time)+x$res2))
data1<-data.frame(data1[[1]])
data1$ID<-as.factor(data1$ID)
data1$group<-as.factor(data1$time)
data1<-data1[order(data1$ID,data1$time),]
head(data1)

d1<-data1[data1$time==-2,]
d2<-data1[data1$time==-1,]
d3<-data1[data1$time==-0,]
d4<-data1[data1$time== 1,]
d5<-data1[data1$time== 2,]

prior1<-list(R=list(V=diag(2)*0.02, nu=3), 
			G=list(G1=list(V=diag(2),
						nu=2,
						alpha.mu=rep(0,2),
						alpha.V=diag(2)*1000)))

m1d<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d1,
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m2d<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d2, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m3d<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d3, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m4d<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d4, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m5d<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d5, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)


## Example 5. two traits, with zero genetic covariance, both with GEI rG 0 to highly positive
# simulate breeding values for 200 genotypes, for two traits. Each trait has intercept variance of 1; one trait has GEI with variance of 0.5; the other has no GEI variance; traits are uncorrelated
G_A<-matrix(c(
   1.0,  -0.3,  -0.9,  -0.64,   
   -0.3,  0.5,  -0.64, -0.45, 
   -0.9,  -0.64,  1.0,  0.3,
   -0.64, -0.45,  0.3,  0.5),4,4)

bv_a <- lapply(1, mvrnorm, n=200, mu=rep(0,4), Sigma=G_A)


# simulate a dataframe from the G matrix 
sire1<-lapply(bv_a, function(x) data.frame(ID=rep(1:nrow(x),each=20), time=rep(-2:2,nrow(x)), intercept1=rep(x[,1],each=20),slope1=rep(x[,2],each=20),
res1=rnorm(4000,0,1.35^2), intercept2=rep(x[,3],each=20),slope2=rep(x[,4],each=20),res2=rnorm(4000,0,1.35^2)))

data1<-lapply(sire1, function(x) cbind(ID=x$ID, time=x$time, phen1=x$intercept1+(x$slope1*x$time)+x$res1, phen2=x$intercept2+(x$slope2*x$time)+x$res2))
data1<-data.frame(data1[[1]])
data1$ID<-as.factor(data1$ID)
data1$group<-as.factor(data1$time)
data1<-data1[order(data1$ID,data1$time),]
head(data1)

d1<-data1[data1$time==-2,]
d2<-data1[data1$time==-1,]
d3<-data1[data1$time==-0,]
d4<-data1[data1$time== 1,]
d5<-data1[data1$time== 2,]

prior1<-list(R=list(V=diag(2)*0.02, nu=3), 
			G=list(G1=list(V=diag(2),
						nu=2,
						alpha.mu=rep(0,2),
						alpha.V=diag(2)*1000)))

m1e<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d1,
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m2e<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d2, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m3e<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d3, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m4e<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d4, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m5e<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d5, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

##Example 6
G_A<-matrix(c(
   1.0,  0.3,  0,  0,   
   0.3,  0.5,  0, 0, 
   0,  0, 1.0,  0.3,
   0, 0, 0.3,  0.5),4,4)

bv_a <- lapply(1, mvrnorm, n=200, mu=rep(0,4), Sigma=G_A)


# simulate a dataframe from the G matrix 
sire1<-lapply(bv_a, function(x) data.frame(ID=rep(1:nrow(x),each=20), time1=rep(-2:2,nrow(x)), time2=rep(2:-2,nrow(x)), intercept1=rep(x[,1],each=20),slope1=rep(x[,2],each=20),
res1=rnorm(4000,0,1.35^2), intercept2=rep(x[,3],each=20),slope2=rep(x[,4],each=20),res2=rnorm(4000,0,1.35^2)))

data1<-lapply(sire1, function(x) cbind(ID=x$ID, time=x$time1, phen1=x$intercept1+(x$slope1*x$time1)+x$res1, phen2=x$intercept2+(x$slope2*x$time2)+x$res2))
data1<-data.frame(data1[[1]])
data1$ID<-as.factor(data1$ID)
data1$group<-as.factor(data1$time)
data1<-data1[order(data1$ID,data1$time),]
head(data1)

d1<-data1[data1$time==-2,]
d2<-data1[data1$time==-1,]
d3<-data1[data1$time==-0,]
d4<-data1[data1$time== 1,]
d5<-data1[data1$time== 2,]

prior1<-list(R=list(V=diag(2)*0.02, nu=3), 
			G=list(G1=list(V=diag(2),
						nu=2,
						alpha.mu=rep(0,2),
						alpha.V=diag(2)*1000)))

m1f<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d1,
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m2f<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d2, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m3f<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d3, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m4f<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d4, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)

m5f<-MCMCglmm(cbind(phen1,phen2) ~ trait-1, 
	random=~us(trait):ID, rcov=~ us(trait):units, data=d5, 
	family=c("gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=50000, burnin=20000, thin=30)



# full cross-env G matrix for each trait
prior2<-list(R=list(V=diag(1)*0.02, nu=0.002), 
			G=list(G1=list(V=diag(5),
						nu=5,
						alpha.mu=rep(0,5),
						alpha.V=diag(5)*1000)))

actualchange_tr1<-MCMCglmm(phen1 ~ group, 
	random=~us(group):ID, data=data1, 
	prior=prior2, verbose=FALSE, nit=15000, burnin=5000, thin=10)

m2<-MCMCglmm(phen2 ~ group, 
	random=~us(group):ID, data=data1, 
	prior=prior1, verbose=FALSE, nit=15000, burnin=5000, thin=10)
	
	
	act_tr1<-data.frame(actualchange_tr1$VCV[,1:(5)^2])
posterior.mode(as.mcmc(act_tr1[,1:25]))




##################################################################################################################################################
##################################################################################################################################################
### 										Analysis King, Roff, Fairburn data																	##
###																																				##
##################################################################################################################################################

data<-read.csv("~/Dropbox/Subspace/RCodes/Working_Matt/Kingetal_data.csv",header=T)
head(data)

data1<-data[data$wng_morph=='L',]
data1<-data.frame(data1$Animal,data1$Sire,data1$Food,data1$DLM_mass,data1$ovary_mass,data1$rem_mass,data1$total_mass,data1$Acquisition, data1$DLMenergy,data1$Ovaryenergy)

names(data1)<-c('Animal','Sire','Food','DLM_mass','ovary_mass','rem_mass','total_mass','Acquisition','DLMenergy','Ovaryenergy')
head(data1)

data1$fa<- (data1$DLMenergy / data1$Acquisition )^2
data1$fa1<-5*((data1$fa - mean(data1$fa)) / sd(data1$fa))

data1$ra<- log(100*(data1$Ovaryenergy / data1$Acquisition )+1)
data1$ra1<-5*((data1$ra - mean(data1$ra)) / sd(data1$ra))

data1$t_off_ac <- data1$DLMenergy + data1$Ovaryenergy
data1$t_off_ac1<- 5*((data1$t_off_ac - mean(data1$t_off_ac)) / sd(data1$t_off_ac))

data1$t_off_al <- (data1$DLMenergy / data1$t_off_ac )^2
data1$t_off_al1<- 5*((data1$t_off_al - mean(data1$t_off_al)) / sd(data1$t_off_al))

data1$tac1<-5*((data1$Acquisition - mean(data1$Acquisition)) / sd(data1$Acquisition))


d1<-data1[data1$Food=='Low',]
d2<-data1[data1$Food=='Med',]
d3<-data1[data1$Food=='AD',]

d1$Sire<-factor(d1$Sire)
d2$Sire<-factor(d2$Sire)
d3$Sire<-factor(d3$Sire)

## re-scale data to enable good model convergence, this will be corrected in later stats
d1$fa1<-d1$fa*10
d1$ra1<-d1$ra*10
d1$t_off_ac1<-d1$t_off_ac*10
d1$t_off_al1<-d1$t_off_al*10
d1$tac1<-d1$tac*10

d2$fa1<-d2$fa*2
d2$ra1<-d2$ra*2
d2$t_off_ac1<-d2$t_off_ac*10
d2$t_off_al1<-d2$t_off_al*10
d2$tac1<-d2$tac*10

d3$fa1<-d3$fa*10
d3$ra1<-d3$ra*10
d3$t_off_ac1<-d3$t_off_ac*10
d3$t_off_al1<-d3$t_off_al*10
d3$tac1<-d3$tac*10



library(MCMCglmm)

prior1<-list(R=list(V=diag(5)*1e-6, nu=6), 
			G=list(G1=list(V=diag(5)*0.5,
						nu=5,
						alpha.mu=rep(0,5),
						alpha.V=diag(5)*100)))


m.low<-MCMCglmm(cbind(fa1, ra1, t_off_ac1, t_off_al1, tac1) ~ trait-1, 
	random=~us(trait):Sire, rcov=~ us(trait):units, data=d1,
	family=c("gaussian","gaussian","gaussian","gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=20000, burnin=10000, thin=10)
	
m.med<-MCMCglmm(cbind(fa1, ra1, t_off_ac1, t_off_al1, tac1) ~ trait-1, 
	random=~us(trait):Sire, rcov=~ us(trait):units, data=d2,
	family=c("gaussian","gaussian","gaussian","gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=220000, burnin=20000, thin=200)

m.high<-MCMCglmm(cbind(fa1, ra1, t_off_ac1, t_off_al1, tac1) ~ trait-1, 
	random=~us(trait):Sire, rcov=~ us(trait):units, data=d3,
	family=c("gaussian","gaussian","gaussian","gaussian","gaussian"), 
	prior=prior1, verbose=TRUE, nit=20000, burnin=10000, thin=10)




##################################################################################################################################################
##################################################################################################################################################
### 										Tensor methods for King, Roff, Fairburn data														##
###																																				##
##################################################################################################################################################
set.seed(123)
# n<-2
# m<-5	
	# a<-data.frame(m1a$VCV[,1:(n)^2])
	# b<-data.frame(m2a$VCV[,1:(n)^2])
	# c<-data.frame(m3a$VCV[,1:(n)^2])
	# d<-data.frame(m4a$VCV[,1:(n)^2])
	# e<-data.frame(m5a$VCV[,1:(n)^2])

n<-5
m<-3
	a<-data.frame(m.low$VCV[,1:(n)^2])
	b<-data.frame(m.med$VCV[,1:(n)^2])
	c<-data.frame(m.high$VCV[,1:(n)^2])

sum_eigvals<-matrix(0,1000,(n*(n+1)/2)+1)
sum_projmat<-list()
sum_eigtenvecs<-list()
sum_eigtenvals<-list()
proj_mat1<-matrix(0,1000,m)
proj_mat2<-matrix(0,1000,m)
eigten_vec1<-matrix(0,1000,n)
eigten_vec2<-matrix(0,1000,n)
eigten_vals1<-matrix(0,1000,n)
eigten_vals2<-matrix(0,1000,n)

for (z in 1:1000){	


########################################################################################################
## 1. Stack the matrices from each MCMC model into one file where ncol = number of traits and nrow = n traits * i environments

	# get variance-covariance posteriors
	m.1<-matrix(as.numeric(a[z,][,1:(n)^2]),n,n)
	m.2<-matrix(as.numeric(b[z,][,1:(n)^2]),n,n)
	m.3<-matrix(as.numeric(c[z,][,1:(n)^2]),n,n)
	# m.4<-matrix(as.numeric(d[z,][,1:(n)^2]),n,n)
	# m.5<-matrix(as.numeric(e[z,][,1:(n)^2]),n,n)


# Gmatrices<-rbind(m.1,m.2,m.3,m.4,m.5)
Gmatrices<-rbind(m.1,m.2,m.3)

########################################################################################################
## 2. Calculate the matrix representation of the covariance tensor (S) 

covcov<-matrix(ncol=n*(n+1)/2,nrow=n*(n+1)/2) # matrix representation of covariance tensor  "S"

# upper left quadrant of S

for (i in 1:n){
  for (j in 1:n){
    covcov[i,j]=cov(Gmatrices[seq(i,n*m,n),i],Gmatrices[seq(j,n*m,n),j])
}
}

# lower left and upper right quadrants of S
for (k in 1:n){
  count=1
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      covcov[k,count+n]=sqrt(2)*cov(Gmatrices[seq(k,n*m,n),k],Gmatrices[seq(i,n*m,n),j])
      covcov[count+n,k]=sqrt(2)*cov(Gmatrices[seq(k,n*m,n),k],Gmatrices[seq(i,n*m,n),j])

count=count+1
}
}
}

# lower right quadrant of S

countx=1
county=1
for (k in 1:(n-1)){
  for (l in (k+1):n){
    for (i in 1:(n-1)){
      for (j in (i+1):n){
        covcov[countx+n,county+n]=2*cov(Gmatrices[seq(k,n*m,n),l],Gmatrices[seq(i,n*m,n),j])
        covcov[county+n,countx+n]=2*cov(Gmatrices[seq(k,n*m,n),l],Gmatrices[seq(i,n*m,n),j])
        county=county+1}}
countx=countx+1
county=1}}

covcov<-as.data.frame(covcov)

########################################################################################################
## 3. Take eigenvectors and eigenvalues of covariance tensor S, and then for each eigenvector of S construct and eigentensor,
## then take eigenvectors and eigenvalues of each eigentensor

eigvecs<-eigen(covcov)$vectors
eigvals<-eigen(covcov)$values
eigten<-matrix(nrow=n*n*(n+1)/2,ncol=n)

for (eval in 1:(n*(n+1)/2)){ #number of eigenvalues of covariance tensor
for (i in 1:n){ #constructing an n*n eigentensor from each eigenvector of S matrix
    eigten[i+n*(eval-1),i]=eigvecs[i,eval] #this loop calculates the diagonal
}

  count=1
  for (i in 1:(n-1)){
   for (j in (i+1):n){  
     eigten[i+n*(eval-1),j]=1/sqrt(2)*eigvecs[count+n,eval]
     eigten[j+n*(eval-1),i]=1/sqrt(2)*eigvecs[count+n,eval]
     count=count+1
}
}
}


# storage matrices of eigentensor eigenvalues and eigenvectors
eigtenvecs<-matrix(nrow=n*n*(n+1)/2, ncol=n)
eigtenvals<-matrix(nrow=n*n*(n+1)/2, ncol=1)

for (i in 1:(n*(n+1)/2)){
  eigtenvecs[((i-1)*n+1):(i*n),]=t(eigen(eigten[((i-1)*n+1):(i*n),])$vectors) # eigenvectors in rows as in Excel
  eigtenvals[((i-1)*n+1):(i*n),1]=eigen(eigten[((i-1)*n+1):(i*n),])$values
}

########################################################################################################\
## 4. Characterise the genetic variance in a single trait combination (y) by projecting y through the G matrix using t(y)%*%G%*%y\

neigten=n*(n+1)/2

projmat<-matrix(ncol=m,nrow=neigten) #matrix containing dot product for each eigentensor with each population delta G
# Cj,k values for population 1 are in the first column, etc.

for (neigs in 1:neigten){
 for (npop in 1:(m)){
  inside=0
  for (i in 1:n){
   for (j in 1:n){
    term=Gmatrices[(npop-1)*n+i,j]*eigten[(neigs-1)*n+i,j]
    inside=inside+term
					}
				}
  		projmat[neigs,npop]=inside
	}
}

sum_eigvals[z,1:(n*(n+1)/2)]<-eigvals
sum_eigvals[z,((n*(n+1)/2)+1)]<-sum(eigvals)
sum_projmat[[z]]<-projmat
sum_eigtenvecs[[z]]<-eigtenvecs
sum_eigtenvals[[z]]<-eigtenvals
proj_mat1[z,]<-projmat[1,]
proj_mat2[z,]<-projmat[2,]
eigten_vec1[z,]<-eigtenvecs[1,]
eigten_vec2[z,]<-eigtenvecs[2,]
eigten_vals1[z,1:n]<-t(eigtenvals[1:n])
eigten_vals2[z,1:n]<-t(eigtenvals[(n+1):(n+n)])

}

library(MCMCglmm)
# 95% CI of eigenvalues of eigentensors and the total of the eigenvalues
HPDinterval(as.mcmc(sum_eigvals[,1:((n*(n+1)/2)+1)]))

HPDinterval(as.mcmc(sum_eigvals[,1]/sum_eigvals[,((n*(n+1)/2)+1)]))
HPDinterval(as.mcmc(sum_eigvals[,2]/sum_eigvals[,((n*(n+1)/2)+1)]))
posterior.mode(as.mcmc(sum_eigvals[,1]/sum_eigvals[,((n*(n+1)/2)+1)]))
posterior.mode(as.mcmc(sum_eigvals[,2]/sum_eigvals[,((n*(n+1)/2)+1)]))


# For the significant eigentensors, the directional change in variance of each trait for each vector
HPDinterval(as.mcmc(eigten_vals1[,1:n]))
HPDinterval(as.mcmc(eigten_vals2[,1:n]))
posterior.mode(as.mcmc(eigten_vals1[,1:n]))
posterior.mode(as.mcmc(eigten_vals2[,1:n]))


# The proportion of variance of each eigentensor explianed by each trait
	pvar1<-matrix(0,1000,n)
		for(i in 1:1000){
				pvar1[i,1]<-eigten_vals1[i,1]^2 / sum(eigten_vals1[i,1:n]^2)
				pvar1[i,2]<-eigten_vals1[i,2]^2 / sum(eigten_vals1[i,1:n]^2)
				pvar1[i,3]<-eigten_vals1[i,3]^2 / sum(eigten_vals1[i,1:n]^2)
				pvar1[i,4]<-eigten_vals1[i,4]^2 / sum(eigten_vals1[i,1:n]^2)
				pvar1[i,5]<-eigten_vals1[i,5]^2 / sum(eigten_vals1[i,1:n]^2)
}
HPDinterval(as.mcmc(pvar1[,1:n]))
posterior.mode(as.mcmc(pvar1[,1:n]))

	pvar2<-matrix(0,1000,n)
		for(i in 1:1000){
				pvar2[i,1]<-eigten_vals2[i,1]^2 / sum(eigten_vals2[i,1:n]^2)
				pvar2[i,2]<-eigten_vals2[i,2]^2 / sum(eigten_vals2[i,1:n]^2)
				pvar2[i,3]<-eigten_vals2[i,3]^2 / sum(eigten_vals2[i,1:n]^2)
				pvar2[i,4]<-eigten_vals2[i,4]^2 / sum(eigten_vals2[i,1:n]^2)
				pvar2[i,5]<-eigten_vals2[i,5]^2 / sum(eigten_vals2[i,1:n]^2)
}
HPDinterval(as.mcmc(pvar2[,1:n]))
posterior.mode(as.mcmc(pvar2[,1:n]))


# Genetic variance of first axis of GEI for each env
HPDinterval(as.mcmc(sqrt(proj_mat1[,1:m]^2)))
posterior.mode(as.mcmc(sqrt(proj_mat1[,1:m]^2)))

# Genetic variance of second axis of GEI for each env...if overlaps zero then no evidnece of any change across env
HPDinterval(as.mcmc(sqrt(proj_mat2[,1:m]^2)))
posterior.mode(as.mcmc(sqrt(proj_mat2[,1:m]^2)))


##################################################################################################################################################
##################################################################################################################################################
### 									Matrix comparison methods for King, Roff, Fairburn data													##
###																																				##
##################################################################################################################################################
library(MCMCglmm)
library(mvtnorm)
library(MASS)
# Gather together the model outputs - VCV components
n<-5
	a<-data.frame(m.low$VCV[,1:(n)^2])
	b<-data.frame(m.med$VCV[,1:(n)^2])
	b<-data.frame(m.high$VCV[,1:(n)^2]) 
	###### NOTE THAT THIS last one IS ALSO LABELLED B? WHY?  Just to make a-b comps below?

low<-matrix(posterior.mode(mcmc(m.low$VCV[,1:(n)^2])),5,5)

h2<-matrix(0,1000,5)
for (i in 1:1000){
	h2mat<-matrix(m.high$VCV[i,1:25],5,5)/10 / ((matrix(m.high$VCV[i,1:25],5,5)/10 + matrix(m.high$VCV[i,26:50],5,5)/10))
	h2[i,]<-diag(h2mat)
}
	
cov2cor(low)
med<-matrix(posterior.mode(mcmc(m.med$VCV[,1:(n)^2])),5,5)
cov2cor(med)
high<-matrix(posterior.mode(mcmc(m.high$VCV[,1:(n)^2])),5,5)
cov2cor(high)


# Collection Zones for each statistic
dist<-matrix(0,1000,1)
distdiff<-matrix(0,1000,1)
vol<-matrix(0,1000,3)
tvar<-matrix(0,1000,3)
pvarGmax<-matrix(0,1000,3)
angle1<-matrix(0,1000,1)
ang1diff<-matrix(0,1000,1)
angle2<-matrix(0,1000,1)
ang2diff<-matrix(0,1000,1)
PeigA<-matrix(0,1000,(n-1))
PeigB<-matrix(0,1000,(n-1))
eigvals1<-matrix(0,1000,(n+1))
eigvals2<-matrix(0,1000,(n+1))
	
# Take two random samples from each model output, 1000 times
a1<-a[sample(nrow(a),1000),] # e1
a2<-a[sample(nrow(a),1000),] # e1
samp1e1<-cbind(a1[,1:(n)^2]) # e1
samp2e1<-cbind(a2[,1:(n)^2]) # e1

a1<-b[sample(nrow(b),1000),] # e2
a2<-b[sample(nrow(b),1000),] # e2
samp1e2<-cbind(a1[,1:(n)^2]) # e2
samp2e2<-cbind(a2[,1:(n)^2]) # e2

# Loop to generate 1000 tests

for (i in 1:1000){
	
	# Create Sampled G-matrices
	d1<-matrix(as.numeric(samp1e1[i,][,1:(n)^2]),n,n)/10
    	d2<-matrix(as.numeric(samp2e1[i,][,1:(n)^2]),n,n)/10
	d3<-matrix(as.numeric(samp1e2[i,][,1:(n)^2]),n,n)/10
    	d4<-matrix(as.numeric(samp2e2[i,][,1:(n)^2]),n,n)/10
	
	# Estimate univariate distributions underlying multivariate space
	dx1<-dmvnorm(d1,rep(0,dim(d1)[1]),d1)
	dx2<-dmvnorm(d2,rep(0,dim(d2)[1]),d2)
	dx3<-dmvnorm(d3,rep(0,dim(d3)[1]),d3)
	dx4<-dmvnorm(d4,rep(0,dim(d4)[1]),d4)
	
	# Oksavainen difference
	e1diff<-mean(sqrt(0.5 * ((dx1 - dx2)^2)/(dx1 + dx2)))
	e2diff<-mean(sqrt(0.5 * ((dx3 - dx4)^2)/(dx3 + dx4)))
	e1e2d<-mean(sqrt(0.5 * ((dx1 - dx3)^2)/(dx1 + dx3))) + mean(sqrt(0.5 * ((dx2 - dx4)^2)/(dx2 + dx4)))
	
	# Oksavainen distance
	dist[i,][1]<-mean(sqrt(0.5 * ((dx1 - dx3)^2)/(dx1 + dx3)))
	distdiff[i,][1]<-(e1diff + e2diff) - e1e2d

	# Eigensystem of matrices
    	env1_1<-eigen(d1, EISPACK=TRUE)
    	env1_2<-eigen(d2, EISPACK=TRUE)
    	env2_1<-eigen(d3, EISPACK=TRUE)
    	env2_2<-eigen(d4, EISPACK=TRUE)

	# Volumes
	
	# Env1
	vol[i,][1]<-((1/2)*pi^2*sqrt(env1_1$values[1])*sqrt(env1_1$values[2])*sqrt(env1_1$values[3])*sqrt(env1_1$values[4])*sqrt(env1_1$values[5]))
	# Env2
 	vol[i,][2]<-((1/2)*pi^2*sqrt(env2_1$values[1])*sqrt(env2_1$values[2])*sqrt(env2_1$values[3])*sqrt(env2_1$values[4])*sqrt(env2_1$values[5]))
	# Env1 - Env2
	vol[i,][3]<-(((1/2)*pi^2*sqrt(env1_1$values[1])*sqrt(env1_1$values[2])*sqrt(env1_1$values[3])*sqrt(env1_1$values[4])*sqrt(env1_1$values[5]))-
				((1/2)*pi^2*sqrt(env2_1$values[1])*sqrt(env2_1$values[2])*sqrt(env2_1$values[3])*sqrt(env2_1$values[4])*sqrt(env2_1$values[5])))
	
	tvar[i,][1]<-(env1_1$values[1]+env1_1$values[2]+env1_1$values[3]+env1_1$values[4]+env1_1$values[5])
	tvar[i,][2]<-(env2_1$values[1]+env2_1$values[2]+env2_1$values[3]+env2_1$values[4]+env2_1$values[5])
	tvar[i,][3]<-(env1_1$values[1]+env1_1$values[2]+env1_1$values[3]+env1_1$values[4]+env1_1$values[5]) - (env2_1$values[1]+env2_1$values[2]+env2_1$values[3]+env2_1$values[4]+env2_1$values[5])
	
	# Variance in Gmax
	pvarGmax[i,][1]<-env1_1$values[1]/sum(env1_1$values)
	pvarGmax[i,][2]<-env2_1$values[1]/sum(env2_1$values)
	pvarGmax[i,][3]<-(env1_1$values[1]/sum(env1_1$values)) - (env2_1$values[1]/sum(env2_1$values))
	
	# Krzanowski test to get angles between Gmax
	e1a1<-krzanowski.test(d1, d2, vecsA=1, vecsB=1)$angles 
	e2a1<-krzanowski.test(d3, d4, vecsA=1, vecsB=1)$angles

	e1a2<-krzanowski.test(d1, d2, vecsA=2, vecsB=2)$angles 
	e2a2<-krzanowski.test(d3, d4, vecsA=2, vecsB=2)$angles

	e1_2a1<-(krzanowski.test(d1, d3, vecsA=1, vecsB=1)$angles + krzanowski.test(d2, d4, vecsA=1, vecsB=1)$angles)	  
	e1_2a2<-(krzanowski.test(d1, d3, vecsA=2, vecsB=2)$angles + krzanowski.test(d2, d4, vecsA=2, vecsB=2)$angles)

	# Angle Stats
	angle1[i,][1]<-(e1_2a1)/2

	ang1diff[i,][1]<-(e1a1 + e2a1) - e1_2a1

	angle2[i,][1]<-(e1_2a2)/2

	ang2diff[i,][1]<-(e1a2 + e2a2) - e1_2a2

		# No Sig Eigenvectors for each matrix
      	
		total1<-sum(env1_1$values)
      	total2<-sum(env1_2$values)
		
		total3<-sum(env2_1$values)
      	total4<-sum(env2_2$values)
		
		for(j in 1:(n-1)){
			idx<-1:j
			tmp1<-env1_1$values[idx]
			tmp2<-env1_2$values[idx]
			
			tmp3<-env1_1$values[idx]
			tmp4<-env1_2$values[idx]
						
			PeigA[i,][j]<-(((sum(tmp1)-sum(tmp2)) + (total1-total2)) - ((sum(tmp1)-total1) + (sum(tmp2) - total2)))
			PeigB[i,][j]<-(((sum(tmp3)-sum(tmp4)) + (total3-total4)) - ((sum(tmp3)-total3) + (sum(tmp4) - total4)))
			}

	eigvals1[i,1:n]<-env1_1$values
	eigvals2[i,1:n]<-env2_1$values

	eigvals1[i,(n+1)]<-sum(env1_1$values)
	eigvals2[i,(n+1)]<-sum(env2_1$values)	
}


# Estimates of the distances between matrices and it's standard deviation
posterior.mode(mcmc(dist[,1]))
HPDinterval(mcmc(dist[,1]))

# Proportion of times out of 1000 where matrices differ in distance
length(distdiff[,1][distdiff[,1]<0])/1000

# Estimates of the difference in volume between the matricies and it's standard deviation
posterior.mode(as.mcmc(vol[,1:3]))
HPDinterval(mcmc(vol[,1:3]))

# Estimates of the difference in total eigenvariance between the matricies and it's standard deviation
posterior.mode(as.mcmc(tvar[,1:3]))
HPDinterval(mcmc(tvar[,1:3]))

# Estimates of the proportion of variance explained by Gmax it's standard deviation
posterior.mode(mcmc(pvarGmax[,1:3]))
HPDinterval(mcmc(pvarGmax[,1:3]))

# Estimates of the angle between the first eigenvector and it's standard deviation
posterior.mode(mcmc(angle1[,1]))
HPDinterval(mcmc(angle1[,1]))

# Proportion of times matrices differ in angle of first eigenvector
1-length(ang1diff[,1][ang1diff[,1]<0])/1000

# Estimates of the angle between the first eigenvector and it's standard deviation
posterior.mode(mcmc(angle2[,1]))
HPDinterval(mcmc(angle2[,1]))

# Proportion of times matrices differ in angle of first eigenvector
1-length(ang2diff[,1][ang2diff[,1]<0])/1000

# Number of significant vectors
1-length(PeigA[,1][PeigA[,1]<0])/1000
1-length(PeigB[,1][PeigB[,1]<0])/1000

# 95% CI of eigenvalues 
HPDinterval(as.mcmc(eigvals1[,1:(n+1)]))
HPDinterval(as.mcmc(eigvals1[,1]/eigvals1[,(n+1)]))
HPDinterval(as.mcmc(eigvals1[,2]/eigvals1[,(n+1)]))
posterior.mode(as.mcmc(eigvals1[,1]/eigvals1[,(n+1)]))
posterior.mode(as.mcmc(eigvals1[,2]/eigvals1[,(n+1)]))


# 95% CI of eigenvalues 
HPDinterval(as.mcmc(eigvals2[,1:(n+1)]))
HPDinterval(as.mcmc(eigvals2[,1]/eigvals2[,(n+1)]))
HPDinterval(as.mcmc(eigvals2[,2]/eigvals2[,(n+1)]))
posterior.mode(as.mcmc(eigvals2[,1]/eigvals2[,(n+1)]))
posterior.mode(as.mcmc(eigvals2[,2]/eigvals2[,(n+1)]))


#Therefore, we have:
#
#1. Number of statistically supported eigenvectors of each matrix: 
#	- informs us as to the level of constraint and number independent traits upon which selection can act.
#	- if fitness measure is included then it informs us as to the strength of directional selection, if many eignevectors then traits unconstrained and unlinked to fitness
#	- the loadings inform us as to the relative strength of selection on each trait - if most variation come from association between one trait and fitness then they will load the highest.
#2. The eigenvalue ratio of the significant eigenvectors describe differences in shape of two matrices
#	- informs us of the amount of variation across ages/environments in each axis, and whether this changes
#3. Angles between significant eigenvectors and statistical support for the number of times the angles differ 
#	- Enables us to compare the orientation of matrices across environments/ages and thus whether the the relationships between traits change
#4. Distance between matrices and statistical support for the number of times the matrices differ
#	- Enables us to test how different the underlying probability distributions are and thus whether there is a different genetic basis across ages/environments, or whether same genes but different effects
#5. Differences in volume between matrices, and the amount that the volumes differ
#	- Informs us as to the amount of additive genetic variation underlying the multivariate distribution, and how that changes across environments.
#
#
#TCI (trait change index): do the traits which define the multi-dimensional space also define those of the other multivariate space
#


