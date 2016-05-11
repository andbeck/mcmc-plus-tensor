#--------------------------------------------------------------------------------------#
## derived.stats function by Robinson and Beckerman 2011 (updated)
## this function takes two MCMCglmm models and the number of traits being modelled
## and returns 12 typical and common comparisons between covariance matrices
##
## the function relies on bayesian posteriors to construct credible intervals
## and or p-value like statistics from the posteriors associated with each comparison
##
## the majority of the tests rely on a comparison of within and between samples from
## the posterior estimates of both matrices (see Ovaskainen 2008)
##
## Now Hosted on GitHub
##
## 11 May 2016 - no longer assumes fixed 1000 samples.  
## Still requires minimum 1000, but more are allowed
## Uses min(m1, m2) of sample size reported by MCMCglmm.summary


#--------------------------------------------------------------------------------------#


derived.stats2<-function(model1,model2,no.traits){
	print("This may take some time..... 30 seconds for 5 traits with ~ 1000 sample size on a 2.4 core2duo macbook")
	if (require(mvtnorm) == FALSE)
        stop("mvtnorm not loaded")
	if (!inherits(model1, "MCMCglmm")) 
        stop("use only with \"MCMCglmm\" objects")
    if (!inherits(model2, "MCMCglmm")) 
        stop("use only with \"MCMCglmm\" objects")
    if (dim(model1$VCV)[1]<1000|dim(model2$VCV)[1]<1000) 
        stop("You need joint posterior dimension > 1000; refit your models")	

	# get variance-covariance posteriors
	a<-data.frame(model1$VCV[,1:(no.traits)^2])
	b<-data.frame(model2$VCV[,1:(no.traits)^2])
	
	# find sample size, set to min if the models have different sample sizes
	rrows_a<-dim(a)[1]
	rrows_b<-dim(b)[1]

	if(rrows_a!=rrows_b){
		rrows=min(rrows_a, rrows_b)
	} else
	{rrows = rrows_a}
	
	# Collection Zones for each statistic
	angle1<-matrix(0,rrows,1)
	ang1diff<-matrix(0,rrows,1)
	angle2<-matrix(0,rrows,1)
	ang2diff<-matrix(0,rrows,1)
	dist<-matrix(0,rrows,1)
	distdiff<-matrix(0,rrows,1)
	sumS<-matrix(0,rrows,1)
	sumSdiff<-matrix(0,rrows,1)
	vol<-matrix(0,rrows,3)
	tvar<-matrix(0,rrows,3)
	ratio<-matrix(0,rrows,3)
	ratiodiff<-matrix(0,rrows,1)
	pvarGmax<-matrix(0,rrows,3)
	tci<-matrix(0,rrows,1)
	PeigA<-matrix(0,rrows,(no.traits-1))
	PeigB<-matrix(0,rrows,(no.traits-1))
	evs<-matrix(NA,rrows,2)
	evols<-matrix(NA,rrows,2)
	
	
	# Take two random samples from each model output, rrows times
	a1<-a[sample(nrow(a),rrows),] # e1
	a2<-a[sample(nrow(a),rrows),] # e1
	samp1e1<-cbind(a1[,1:(no.traits)^2]) # e1
	samp2e1<-cbind(a2[,1:(no.traits)^2]) # e1
	
	a1<-b[sample(nrow(b),rrows),] # e2
	a2<-b[sample(nrow(b),rrows),] # e2
	samp1e2<-cbind(a1[,1:(no.traits)^2]) # e2
	samp2e2<-cbind(a2[,1:(no.traits)^2]) # e2

	#--------------------------------------------------------------------------------------#
	# Loop to generate rrows > 1000 tests
	#--------------------------------------------------------------------------------------#

	for (i in 1:rrows){
		
		# Create Sampled G-matrices
		d1<-matrix(as.numeric(samp1e1[i,][,1:(no.traits)^2]),no.traits,no.traits)
	    	d2<-matrix(as.numeric(samp2e1[i,][,1:(no.traits)^2]),no.traits,no.traits)
		d3<-matrix(as.numeric(samp1e2[i,][,1:(no.traits)^2]),no.traits,no.traits)
	    	d4<-matrix(as.numeric(samp2e2[i,][,1:(no.traits)^2]),no.traits,no.traits)
		
		# Estimate univariate distributions underlying multivariate space
		dx1<-dmvnorm(d1,rep(0,dim(d1)[1]),d1)
		dx2<-dmvnorm(d2,rep(0,dim(d2)[1]),d2)
		dx3<-dmvnorm(d3,rep(0,dim(d3)[1]),d3)
		dx4<-dmvnorm(d4,rep(0,dim(d4)[1]),d4)
		
		# Oksavainen differences
		e1diff<-mean(sqrt(0.5 * ((dx1 - dx2)^2)/(dx1 + dx2)))
		e2diff<-mean(sqrt(0.5 * ((dx3 - dx4)^2)/(dx3 + dx4)))
		e1e2d<-mean(sqrt(0.5 * ((dx1 - dx3)^2)/(dx1 + dx3))) + mean(sqrt(0.5 * ((dx2 - dx4)^2)/(dx2 + dx4)))
		
		# Oksavainen distance test
		dist[i,][1]<-mean(sqrt(0.5 * ((dx1 - dx3)^2)/(dx1 + dx3)))
		distdiff[i,][1]<-(e1diff + e2diff) - e1e2d
	
		# Krzanowski tests to get angles between Gmax
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
	
		# Eigensystem of matrices 
	    	env1_1<-eigen(d1)
	    	env1_2<-eigen(d2)
	    	env2_1<-eigen(d3)
	    	env2_2<-eigen(d4)
	
		# No Sig Eigenvectors for each matrix
      	
		total1<-sum(env1_1$values)
      	total2<-sum(env1_2$values)
		
		total3<-sum(env2_1$values)
      	total4<-sum(env2_2$values)
		
		for(j in 1:(no.traits-1)){
			idx<-1:j
			tmp1<-env1_1$values[idx]
			tmp2<-env1_2$values[idx]
			
			tmp3<-env2_1$values[idx]
			tmp4<-env2_2$values[idx]
						
			PeigA[i,][j]<-(((sum(tmp1)-sum(tmp2)) + (total1-total2)) - ((sum(tmp1)-total1) + (sum(tmp2) - total2)))
			PeigB[i,][j]<-(((sum(tmp3)-sum(tmp4)) + (total3-total4)) - ((sum(tmp3)-total3) + (sum(tmp4) - total4)))
			}
		
		### Volumes
		# Two types
		# prob dist style (pi^2*sqrtvals)
		# sum eigenvalues
		#####
		
		#prob dist style
		sqrt_vvalues1<-sqrt(env1_1$values)
		sqrt_vvalues2<-sqrt(env2_1$values)
		
		vol1<-(1/2)*pi^2*prod(sqrt_vvalues1)
		vol2<-(1/2)*pi^2*prod(sqrt_vvalues2)
		
		vol[i,][1]<-vol1 # allocate E1
		vol[i,][2]<-vol2 # allocate E2
	 	vol[i,][3]<-vol1-vol2 # allocate difference
		
		# sum style
		tvar[i,][1]<-sum(env1_1$values) # allocate E1
		tvar[i,][2]<-sum(env2_1$values) # allocate E2
		tvar[i,][3]<-sum(env1_1$values) - sum(env2_1$values) # allocate differences
		
		# First and Second EV Ratios: Eccentricity
		ratio[i,][1]<-(env1_1$values[1]/env1_1$values[2])
		ratio[i,][2]<-(env2_1$values[1]/env2_1$values[2])
		ratio[i,][3]<-((env1_1$values[1]/env1_1$values[2])-(env2_1$values[1]/env2_1$values[2]))
	
		#Eccentricity diff
		ratiodiff[i,][1]<- ((((env1_1$values[1]/env1_1$values[2]) - 
							(env1_2$values[1]/env1_2$values[2])) + 
							((env2_1$values[1]/env2_1$values[2]) - 
							(env2_2$values[1]/env2_2$values[2])))-
							(((env1_1$values[1]/env1_1$values[2]) - 
							(env2_1$values[1]/env2_1$values[2])) + 
							((env1_2$values[1]/env1_2$values[2]) - 
							(env2_2$values[1]/env2_2$values[2]))))
	
		# Variance in Gmax
		pvarGmax[i,][1]<-env1_1$values[1]/sum(env1_1$values)
		pvarGmax[i,][2]<-env2_1$values[1]/sum(env2_1$values)
		pvarGmax[i,][3]<-env1_1$values[1]/sum(env1_1$values) - env2_1$values[1]/sum(env2_1$values)
		
		# Kirkpatrick Stats
		
		# ND
		evs1<-sum(env1_1$values)/env1_1$values[1]
		evs2<-sum(env2_1$values)/env2_1$values[1]

		# Evolvability
		evol1<-sqrt(env1_1$values[1])
		evol2<-sqrt(env2_1$values[1])

		evs[i,]<-c(evs1,evs2)
		evols[i,]<-c(evol1,evol2)
		
		# Trait Change Index
	
		# stuff needed	
		E1Vals<-env1_1$values
		E1Vecs<-env1_1$vectors[,1:3]
		E2Vals<-env2_1$values
	
		# calculate prop in threespace for E1
		E1_map<-round(100*sum(E1Vals[1:3])/sum(E1Vals))
		# project E2 into E1 space
		E2proj<-t(E1Vecs) %*% d3 %*% (E1Vecs)
		# calc eigensystem of projected E2
		e2Proj.eig<-eigen(E2proj)
		E2ProjVals<-e2Proj.eig$values
		#calculate prop in threespace for E2
		E2proj_map<-round(100*sum(E2ProjVals[1:3])/sum(E2Vals))
		
		# calcuate difference in proportion - if includes 0, no change in traits contributing to G.
		tci[i,][1] <-E1_map - E2proj_map	
		
	}
	
	#--------------------------------------------------------------------------------------#
	# collection and printing of all stats
	# (a 3 column matrix with labels)
	#--------------------------------------------------------------------------------------#

	#--------------------------------------------------------------------------------------#
	## vec calc to help presentation
	#--------------------------------------------------------------------------------------#
	A<-sum(apply(PeigA,2,function(x)length(x[x<0]))<50)
	B<-sum(apply(PeigB,2,function(x)length(x[x<0]))<50)
	
	#--------------------------------------------------------------------------------------#
	## collection and printing
	#--------------------------------------------------------------------------------------#
	mcore<-matrix(NA,12,4,
		dimnames=list(c(
			"Ovaskainen D",
			"Variance Gmax 1",
			"Variance Gmax 2",
			"VarGmax Diff",
			"Angle Between Gmax",
			"prob-Volume 1",
			"prob-Volume 2",
			"prob-VolDiff",
			"sum-Volume 1",
			"sum-Volume 2",
			"sum-VolDiff",
			"TCI"),
			c("Mode/Difference", "Lower","Upper", "P[psi<0|p<0.05]*")))		

	mvecs<-matrix(NA,1,2,
		dimnames=list(c("No. Vectors"),	c("Matrix A", "Matrix B")))		
		
	Nd_evol<-matrix(NA,4,3,
			dimnames=list(c("NdA","NdB","EvolA","EvolB"),c("Mode/Difference", "Lower","Upper")))
	
	meccs<-matrix(NA,3,4,
			dimnames=list(c("Eccentricity 1","Eccentricity 2","EccDiff"),
			c("Mode/Difference", "Lower","Upper","P[psi<0|p<0.05]*")))
	
	
	## CORE	
	#Ov D	
	mcore[1,1:4]<-c(posterior.mode(mcmc(dist[,1])), HPDinterval(mcmc(dist[,1])),(1-length(distdiff[,1][distdiff[,1]<0])/rrows))
	
	# Gmax		
	mcore[2,1:4]<-c(posterior.mode(mcmc(pvarGmax[,1])), HPDinterval(mcmc(pvarGmax[,1])),
		ifelse(all(HPDinterval(mcmc(pvarGmax[,1]))>0)|all(HPDinterval(mcmc(pvarGmax[,1]))<0),0.05,NA))		# Variance in Gmax (1st axis) in E1
	mcore[3,1:4]<-c(posterior.mode(mcmc(pvarGmax[,2])), HPDinterval(mcmc(pvarGmax[,2])),
		ifelse(all(HPDinterval(mcmc(pvarGmax[,2]))>0)|all(HPDinterval(mcmc(pvarGmax[,2]))<0),0.05,NA))		# Variance in Gmax (1st axis) in E2
	mcore[4,1:4]<-c(posterior.mode(mcmc(pvarGmax[,3])), HPDinterval(mcmc(pvarGmax[,3])),
		ifelse(all(HPDinterval(mcmc(pvarGmax[,3]))>0)|all(HPDinterval(mcmc(pvarGmax[,3]))<0),0.05,NA))

	# Angles		
	mcore[5,1:4]<-c(posterior.mode(mcmc(angle1[,1])), HPDinterval(mcmc(angle1[,1])), (1-length(ang1diff[,1][ang1diff[,1]<0])/rrows))
			
	# Difference in prob-Volume (if all positive, A>B; if all negative B>A)
	mcore[6,1:4]<-c(posterior.mode(mcmc(vol[,1])),HPDinterval(mcmc(vol[,1])),
		ifelse(all(HPDinterval(mcmc(vol[,1]))>0)|all(HPDinterval(mcmc(vol[,1]))<0),0.05,NA))					
	mcore[7,1:4]<-c(posterior.mode(mcmc(vol[,2])),HPDinterval(mcmc(vol[,2])),
		ifelse(all(HPDinterval(mcmc(vol[,2]))>0)|all(HPDinterval(mcmc(vol[,2]))<0),0.05,NA))					
	mcore[8,1:4]<-c(posterior.mode(mcmc(vol[,3])),HPDinterval(mcmc(vol[,3])),
		ifelse(all(HPDinterval(mcmc(vol[,3]))>0)|all(HPDinterval(mcmc(vol[,3]))<0),0.05,NA))

	# Difference in sum-Volume (if all positive, A>B; if all negative B>A)
	mcore[9,1:4]<-c(posterior.mode(mcmc(tvar[,1])),HPDinterval(mcmc(tvar[,1])),
		ifelse(all(HPDinterval(mcmc(tvar[,1]))>0)|all(HPDinterval(mcmc(tvar[,1]))<0),0.05,NA))					
	mcore[10,1:4]<-c(posterior.mode(mcmc(tvar[,2])),HPDinterval(mcmc(tvar[,2])),
		ifelse(all(HPDinterval(mcmc(tvar[,2]))>0)|all(HPDinterval(mcmc(tvar[,2]))<0),0.05,NA))					
	mcore[11,1:4]<-c(posterior.mode(mcmc(tvar[,3])),HPDinterval(mcmc(tvar[,3])),
		ifelse(all(HPDinterval(mcmc(tvar[,3]))>0)|all(HPDinterval(mcmc(tvar[,3]))<0),0.05,NA))
		
	# Trait Change Index - Proportion of Variance in matrix 2 exp by 1
	mcore[12,1:4]<-c(posterior.mode(mcmc(tci[,1])),HPDinterval(mcmc(tci[,1])),
		ifelse(all(HPDinterval(mcmc(tci[,1]))>0)|all(HPDinterval(mcmc(tci[,1]))<0),0.05,NA))				
	
	
	## VECTORS
	mvecs[1,]<-c(ifelse(A==0,1,A),ifelse(B==0,1,B))

	## Kirkpatrick
	Nd_evol[1,]<-c(posterior.mode(mcmc(evs[,1])),HPDinterval(mcmc(evs[,1])))
		Nd_evol[2,]<-c(posterior.mode(mcmc(evs[,2])),HPDinterval(mcmc(evs[,2])))
			Nd_evol[3,]<-c(posterior.mode(mcmc(evols[,1])),HPDinterval(mcmc(evols[,1])))
				Nd_evol[4,]<-c(posterior.mode(mcmc(evols[,2])),HPDinterval(mcmc(evols[,2])))
	
	## ECCENTRICITY	(Ratio of first to second EV)
	meccs[1,1:4]<-c(posterior.mode(mcmc(ratio[,1])),HPDinterval(mcmc(ratio[,1])),
				ifelse(all(HPDinterval(mcmc(ratio[,1]))>0)|all(HPDinterval(mcmc(ratio[,1]))<0),0.05,NA))
	meccs[2,1:4]<-c(posterior.mode(mcmc(ratio[,2])),HPDinterval(mcmc(ratio[,2])),
				ifelse(all(HPDinterval(mcmc(ratio[,2]))>0)|all(HPDinterval(mcmc(ratio[,2]))<0),0.05,NA))
	meccs[3,1:4]<-c(posterior.mode(mcmc(ratio[,3])),HPDinterval(mcmc(ratio[,3])), (1-length(ratiodiff[,1][ratiodiff[,1]<0])/rrows))						
	
	cat(paste("Here are the derived Statistics - some acos() warnings are expected out of ",rrows, "iterations",
		"\n","TCI will be 0 if there are only 3 traits","\n",
		"*[psi<0] is appropriate for Ovaskainen D, Angle and Vectors. All others are 95% Credible Intervals and p<0.05 indicates that the CI does not include 0.","\n","\n","*If psi<0.05,between differences were larger than within differences 95% of ", rrows, "samples","\n","\n"))
	return(list(CoreStats=round(mcore,3),Vectors=mvecs,Nd_Evolv=round(Nd_evol,3),Eccentricity=round(meccs,3)))

}
