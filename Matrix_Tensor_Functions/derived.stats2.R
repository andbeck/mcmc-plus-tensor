# --------------------------------------------------------------------------------------#
## derived.stats function by Robinson and Beckerman 2011 Ecology Letters (updated)
##
## this function takes two MCMCglmm models and the number of traits being modelled
## and returns 12 typical and common comparisons between covariance matrices
## output is in the form of a list with three slots
## Core Stats (D, angles, volumes), Vector Stats (vectors, kirkpatrick stats) and Eccentricity Stats
##
## the function relies on bayesian posteriors to construct credible intervals
## and or p-value like statistics from the posteriors associated with each comparison
##
## the D-divergence (Ovaskainen) and angles (krzanowski) tests rely on a comparison of 
## within and between samples from the posterior estimates of both matrices 
## as defined in (see Ovaskainen 2008)
##
## UPDATED 11 May 2016 - no longer assumes fixed 1000 samples.  
## Still requires minimum 1000, but more are allowed
## Uses min(m1, m2) of sample size reported by MCMCglmm.summary

## UPDATED 4 July 2018 - Ovaskainen's Divergence test was incorrect.  Now fixed.
## Added use of dplyr code in several places.


#--------------------------------------------------------------------------------------#


derived.stats2 <- function(model1,model2,no.traits){
  print("This takes some time..... 30 seconds for 5 traits with ~ 1000 sample size on a 2.4 core2duo macbook")
  
  # package dependencies
  if (require(mvtnorm) == FALSE)
    stop("mvtnorm required and not loaded")
  if (require(dplyr) == FALSE)
    stop("dplyr required and not loaded")
  if(require(MCMCglmm) == FALSE)
    stop("MCMCglmm required and not loaded")
  
  # model dependencies
  if (!inherits(model1, "MCMCglmm")) 
    stop("use only with \"MCMCglmm\" objects")
  if (!inherits(model2, "MCMCglmm")) 
    stop("use only with \"MCMCglmm\" objects")
  
  # posterior dependencies
  if (dim(model1$VCV)[1]<1000|dim(model2$VCV)[1]<1000) 
    stop("You probably want joint posterior dimension > 1000; refit your models")	
  
  # get variance-covariance posteriors
  a <- data.frame(model1$VCV[,1:(no.traits)^2])
  b <- data.frame(model2$VCV[,1:(no.traits)^2])
  
  # find posterior distribution sample size, set to min if the models have different sample sizes
  rrows_a <- dim(a)[1]
  rrows_b <- dim(b)[1]
  
  if(rrows_a != rrows_b){
    rrows = min(rrows_a, rrows_b)
  } else
  {rrows = rrows_a}
  
  # Collection Zones for each statistic
  angle1 <- matrix(0,rrows,1)
  ang1diff <- matrix(0,rrows,1)
  angle2 <- matrix(0,rrows,1)
  ang2diff <- matrix(0,rrows,1)
  
  dist <- matrix(0,rrows,1)
  distdiff <- matrix(0,rrows,1)
  
  sumS <- matrix(0,rrows,1)
  sumSdiff <- matrix(0,rrows,1)
  vol <- matrix(0,rrows,3)
  
  tvar <- matrix(0,rrows,3)
  
  ratio <- matrix(0,rrows,3)
  ratiodiff <- matrix(0,rrows,1)
  
  pvarGmax <- matrix(0,rrows,3)
  
  tci <- matrix(0,rrows,1)
  
  PeigA <- matrix(0,rrows,(no.traits-1))
  PeigB <- matrix(0,rrows,(no.traits-1))
  
  evs <- matrix(NA,rrows,2)
  evols <- matrix(NA,rrows,2)
  
  # scrambled/sampled VCVs: 2 within each model/environment
  sample_1_mod_e1  <-  sample_n(a, rrows) # model 1 sample 1
  sample_2_mod_e1  <-  sample_n(a, rrows) # model 1 sample 2
  sample_1_mod_e2  <-  sample_n(b, rrows) # model 1 sample 1
  sample_2_mod_e2  <-  sample_n(b, rrows) # model 1 sample 2
  
  #--------------------------------------------------------------------------------------#
  # Loop to generate rrows >= 1000 tests
  #--------------------------------------------------------------------------------------#
  
  for (i in 1:rrows){
    
    # build matrices of each random sample of the posteriors ----
    
    # model 1 - sample 1
    df1_mat_1  <-  matrix(unlist(sample_1_mod_e1[i,]), no.traits, no.traits)
    # model 1 - sample 2
    df1_mat_2  <-  matrix(unlist(sample_2_mod_e1[i,]), no.traits, no.traits)
    # model 2 - sample 1
    df2_mat_1  <-  matrix(unlist(sample_1_mod_e2[i,]), no.traits, no.traits)
    # model 2 - sample 2
    df2_mat_2  <-  matrix(unlist(sample_2_mod_e2[i,]), no.traits, no.traits)
    
    # Ovaskainen's D divergence test -------------------------------
    
    # uses Ddivergence() function from MCMCglmm version > 2.6
    # this was updated in MCMCglmm in July 2018
    
    # 2 x within divergences
    # same models, separate samples
    AA  <-  Ddivergence(df1_mat_1, df1_mat_2, n = 1000)
    BB  <-  Ddivergence(df2_mat_1, df2_mat_2, n = 1000)
    
    # 2 x between divergences
    # different models
    AB  <-  Ddivergence(df1_mat_1, df2_mat_1, n = 1000)
    BA  <-  Ddivergence(df1_mat_2, df2_mat_2, n = 1000)
    
    # Divergence Stat
    # could be average of the two (we do this with the angles)
    dist[i,]  <-  AB # divergence statistic (choose one between divergence)
    #dist[i,]  <-  (AB + BA)/2 # divergence statistic (average between divergence)
    
    # Divergence test stat (psi-stat)
    distdiff[i,]  <-  (AA + BB) - (AB + BA) # divergence test statistics (Ovaskainen's psi)
    
    # Krzanowski tests to get angles between Gmax --------------------------
    
    # vector 1 angles for model_1 samples and model 2 samples (environment 1 and 2)
    # these are the within environment (within model) samples
    e1a1 <- krzanowski.test(df1_mat_1, df1_mat_2, vecsA=1, vecsB=1)$angles 
    e2a1 <- krzanowski.test(df2_mat_1, df2_mat_2, vecsA=1, vecsB=1)$angles
    
    # vector 2 angles for model_1 and model_2 (environment 1 and 2)
    # these are the within environment (within model) samples
    e1a2 <- krzanowski.test(df1_mat_1, df1_mat_2, vecsA=2, vecsB=2)$angles 
    e2a2 <- krzanowski.test(df2_mat_1, df2_mat_2, vecsA=2, vecsB=2)$angles
    
    # sum of vector 1 angles for model_1 versus model_2 (2 samples each)
    # these are the between environment (between model) samples
    e1_2a1 <- (krzanowski.test(df1_mat_1, df2_mat_1, vecsA=1, vecsB=1)$angles + 
               krzanowski.test(df1_mat_2, df2_mat_2, vecsA=1, vecsB=1)$angles)	  
    
    # sum of vector 2 angles for model_1 versus model_2 (2 samples each)
    # these are the between environment (between model) samples
    e1_2a2 <- (krzanowski.test(df1_mat_1, df2_mat_1, vecsA=2, vecsB=2)$angles + 
               krzanowski.test(df1_mat_2, df2_mat_2, vecsA=2, vecsB=2)$angles)	  
    
    # Angle 1
    angle1[i,][1] <- (e1_2a1)/2 # average angle between environments
    # Angle 1 test statistic defined by thinking about how angles work and circles
    ang1diff[i,][1] <- (e1a1 + e2a1) - e1_2a1 
    
    # Angle 2
    angle2[i,][1] <- (e1_2a2)/2
    # Angle 2 test statistic
    ang2diff[i,][1] <- (e1a2 + e2a2) - e1_2a2
    
    # No Sig Eigenvectors for each matrix ------------------------
    
    # define eigensystem of matrices
    
    env1_1 <- eigen(df1_mat_1)
    env1_2 <- eigen(df1_mat_2)
    env2_1 <- eigen(df2_mat_1)
    env2_2 <- eigen(df2_mat_2)
    
    # calculate the total variance in matrix  	
    total1 <- sum(env1_1$values)
    total2 <- sum(env1_2$values)
    
    total3 <- sum(env2_1$values)
    total4 <- sum(env2_2$values)
    
    # loop over number of vectors (max = number of traits) 
    # and test for when allowing more vectors 
    # fails to deviate from total variance
    
    for(j in 1:(no.traits-1)){
      idx <- 1:j
      tmp1 <- env1_1$values[idx]
      tmp2 <- env1_2$values[idx]
      
      tmp3 <- env2_1$values[idx]
      tmp4 <- env2_2$values[idx]
      
      PeigA[i,][j] <- (((sum(tmp1)-sum(tmp2)) + (total1-total2)) - ((sum(tmp1)-total1) + (sum(tmp2) - total2)))
      PeigB[i,][j] <- (((sum(tmp3)-sum(tmp4)) + (total3-total4)) - ((sum(tmp3)-total3) + (sum(tmp4) - total4)))
    }
    
    ### Volume statistics ----------------------------------------------
    
    # Two types
    # prob dist style (pi^2*sqrtvals)
    # sum eigenvalues
    # statistics here are simply differences and credible intervals
    
    #prob dist method ----
    sqrt_vvalues1 <- sqrt(env1_1$values)
    sqrt_vvalues2 <- sqrt(env2_1$values)
    
    vol1 <- (1/2)*pi^2*prod(sqrt_vvalues1)
    vol2 <- (1/2)*pi^2*prod(sqrt_vvalues2)
    
    vol[i,][1] <- vol1 # allocate E1
    vol[i,][2] <- vol2 # allocate E2
    vol[i,][3] <- vol1 - vol2 # allocate difference
    
    # sum method ----
    tvar[i,][1] <- sum(env1_1$values) # allocate E1
    tvar[i,][2] <- sum(env2_1$values) # allocate E2
    tvar[i,][3] <- sum(env1_1$values) - sum(env2_1$values) # allocate differences
    
    # First and Second EV Ratios: Eccentricity ---------------------------------------------------
    ratio[i,][1] <- (env1_1$values[1]/env1_1$values[2])
    ratio[i,][2] <- (env2_1$values[1]/env2_1$values[2])
    ratio[i,][3] <- ((env1_1$values[1]/env1_1$values[2])-(env2_1$values[1]/env2_1$values[2]))
    
    #Eccentricity diff
    ratiodiff[i,][1] <-  ((((env1_1$values[1]/env1_1$values[2]) - 
                            (env1_2$values[1]/env1_2$values[2])) + 
                           ((env2_1$values[1]/env2_1$values[2]) - 
                              (env2_2$values[1]/env2_2$values[2])))-
                          (((env1_1$values[1]/env1_1$values[2]) - 
                              (env2_1$values[1]/env2_1$values[2])) + 
                             ((env1_2$values[1]/env1_2$values[2]) - 
                                (env2_2$values[1]/env2_2$values[2]))))
    
    # Variance in Gmax ---------------------------------------------------
    
    pvarGmax[i,][1] <- env1_1$values[1]/sum(env1_1$values)
    pvarGmax[i,][2] <- env2_1$values[1]/sum(env2_1$values)
    pvarGmax[i,][3] <- env1_1$values[1]/sum(env1_1$values) - env2_1$values[1]/sum(env2_1$values)
    
    # Kirkpatrick Stats ---------------------------------------------------
    
    # ND
    evs1 <- sum(env1_1$values)/env1_1$values[1]
    evs2 <- sum(env2_1$values)/env2_1$values[1]
    
    # Evolvability
    evol1 <- sqrt(env1_1$values[1])
    evol2 <- sqrt(env2_1$values[1])
    
    evs[i,] <- c(evs1,evs2)
    evols[i,] <- c(evol1,evol2)
    
    # Trait Change Index ---------------------------------------------------
    # uses projection technique to ask how much variation the trait mapping in E1 explains 
    # in E2. 
    # stuff needed	
    E1Vals <- env1_1$values
    E1Vecs <- env1_1$vectors[,1:3]
    E2Vals <- env2_1$values
    
    # calculate prop in threespace for E1
    E1_map <- round(100*sum(E1Vals[1:3])/sum(E1Vals))
    
    # project E2 into E1 space
    E2proj <- t(E1Vecs) %*% df2_mat_1 %*% (E1Vecs) # df2_mat_1 was d3
    
    # calc eigensystem of projected E2
    e2Proj.eig <- eigen(E2proj)
    E2ProjVals <- e2Proj.eig$values
    
    #calculate prop in threespace for E2
    E2proj_map <- round(100*sum(E2ProjVals[1:3])/sum(E2Vals))
    
    # calcuate difference in proportion - if includes 0, no change in traits 
    # contributing to G.
    tci[i,][1]  <- E1_map - E2proj_map	
    
  }
  
  #--------------------------------------------------------------------------------------#
  # collection and printing of all stats
  # (a 3 column matrix with labels)
  #--------------------------------------------------------------------------------------#
  
  #--------------------------------------------------------------------------------------#
  ## vec calc to help presentation
  #--------------------------------------------------------------------------------------#
  A <- sum(apply(PeigA, 2, function(x) length(x[x<0])) < 50)
  B <- sum(apply(PeigB, 2, function(x) length(x[x<0])) < 50)
  
  #--------------------------------------------------------------------------------------#
  ## collection zone setup ----------
  #--------------------------------------------------------------------------------------#
  mcore <- matrix(NA,12,4,
                dimnames=list(c(
                  "Ovaskainen D ->",
                  "Variance Gmax 1",
                  "Variance Gmax 2",
                  "VarGmax Diff ->",
                  "Angle Between Gmax ->",
                  "prob-Volume 1",
                  "prob-Volume 2",
                  "prob-VolDiff ->",
                  "sum-Volume 1",
                  "sum-Volume 2",
                  "sum-VolDiff->",
                  "TCI ->"),
                  c("Mode/Difference", "Lower","Upper", "P[psi<0|p<0.05]*")))		
  
  mvecs <- matrix(NA,1,2,
                dimnames=list(c("No. Vectors"),	c("Matrix A", "Matrix B")))		
  
  Nd_evol <- matrix(NA,4,3,
                  dimnames=list(c("NdA","NdB","EvolA","EvolB"),c("Mode/Difference", "Lower","Upper")))
  
  meccs <- matrix(NA,3,4,
                dimnames=list(c("Eccentricity 1","Eccentricity 2","EccDiff"),
                              c("Mode/Difference", "Lower","Upper","P[psi<0|p<0.05]*")))
  
  
  ## CORE	Statistics to Present ---------
  #Ov D	
  mcore[1,1:4] <- c(posterior.mode(mcmc(dist[,1])), # estimated D stat
                  HPDinterval(mcmc(dist[,1])), # interval on D stat
                  (1-sum(distdiff<0)/rrows)) # prop within vs between > 0; needs to be smaller than 0.05 for >95% to be less than 0; changed 2018
  
  # Gmax
  
  # Variance in Gmax (1st axis) in E1
  mcore[2,1:4] <- c(posterior.mode(mcmc(pvarGmax[,1])), 
                  HPDinterval(mcmc(pvarGmax[,1])),
                  ifelse(all(HPDinterval(mcmc(pvarGmax[,1]))>0) | all(HPDinterval(mcmc(pvarGmax[,1]))<0), 0.05, NA))		
  
  # Variance in Gmax (1st axis) in E2
  mcore[3,1:4] <- c(posterior.mode(mcmc(pvarGmax[,2])),
                  HPDinterval(mcmc(pvarGmax[,2])),
                  ifelse(all(HPDinterval(mcmc(pvarGmax[,2]))>0) | all(HPDinterval(mcmc(pvarGmax[,2]))<0), 0.05, NA))		
  # Gmax difference in Variance
  mcore[4,1:4] <- c(posterior.mode(mcmc(pvarGmax[,3])), 
                  HPDinterval(mcmc(pvarGmax[,3])),
                  ifelse(all(HPDinterval(mcmc(pvarGmax[,3]))>0) | all(HPDinterval(mcmc(pvarGmax[,3]))<0), 0.05, NA))
  
  # Angles		
  mcore[5,1:4] <- c(posterior.mode(mcmc(angle1[,1])), HPDinterval(mcmc(angle1[,1])), 
                  (1-length(ang1diff[,1][ang1diff[,1]<0])/rrows))
  
  # Difference in prob-Volume (if all positive, A>B; if all negative B>A)
  mcore[6,1:4] <- c(posterior.mode(mcmc(vol[,1])),HPDinterval(mcmc(vol[,1])),
                  ifelse(all(HPDinterval(mcmc(vol[,1]))>0)|all(HPDinterval(mcmc(vol[,1]))<0),0.05,NA))					
  mcore[7,1:4] <- c(posterior.mode(mcmc(vol[,2])),HPDinterval(mcmc(vol[,2])),
                  ifelse(all(HPDinterval(mcmc(vol[,2]))>0)|all(HPDinterval(mcmc(vol[,2]))<0),0.05,NA))					
  mcore[8,1:4] <- c(posterior.mode(mcmc(vol[,3])),HPDinterval(mcmc(vol[,3])),
                  ifelse(all(HPDinterval(mcmc(vol[,3]))>0)|all(HPDinterval(mcmc(vol[,3]))<0),0.05,NA))
  
  # Difference in sum-Volume (if all positive, A>B; if all negative B>A)
  mcore[9,1:4] <- c(posterior.mode(mcmc(tvar[,1])),HPDinterval(mcmc(tvar[,1])),
                  ifelse(all(HPDinterval(mcmc(tvar[,1]))>0)|all(HPDinterval(mcmc(tvar[,1]))<0),0.05,NA))					
  mcore[10,1:4] <- c(posterior.mode(mcmc(tvar[,2])),HPDinterval(mcmc(tvar[,2])),
                   ifelse(all(HPDinterval(mcmc(tvar[,2]))>0)|all(HPDinterval(mcmc(tvar[,2]))<0),0.05,NA))					
  mcore[11,1:4] <- c(posterior.mode(mcmc(tvar[,3])),HPDinterval(mcmc(tvar[,3])),
                   ifelse(all(HPDinterval(mcmc(tvar[,3]))>0)|all(HPDinterval(mcmc(tvar[,3]))<0),0.05,NA))
  
  # Trait Change Index - Proportion of Variance in matrix 2 exp by 1
  mcore[12,1:4] <- c(posterior.mode(mcmc(tci[,1])),HPDinterval(mcmc(tci[,1])),
                   ifelse(all(HPDinterval(mcmc(tci[,1]))>0)|all(HPDinterval(mcmc(tci[,1]))<0),0.05,NA))				
  
  
  ## VECTORS
  mvecs[1,] <- c(ifelse(A == 0, 1, A), ifelse(B == 0, 1, B))
  
  ## Kirkpatrick
  Nd_evol[1,] <- c(posterior.mode(mcmc(evs[,1])),HPDinterval(mcmc(evs[,1])))
  Nd_evol[2,] <- c(posterior.mode(mcmc(evs[,2])),HPDinterval(mcmc(evs[,2])))
  Nd_evol[3,] <- c(posterior.mode(mcmc(evols[,1])),HPDinterval(mcmc(evols[,1])))
  Nd_evol[4,] <- c(posterior.mode(mcmc(evols[,2])),HPDinterval(mcmc(evols[,2])))
  
  ## ECCENTRICITY	(Ratio of first to second EV)
  meccs[1,1:4] <- c(posterior.mode(mcmc(ratio[,1])),HPDinterval(mcmc(ratio[,1])),
                  ifelse(all(HPDinterval(mcmc(ratio[,1]))>0)|all(HPDinterval(mcmc(ratio[,1]))<0),0.05,NA))
  meccs[2,1:4] <- c(posterior.mode(mcmc(ratio[,2])),HPDinterval(mcmc(ratio[,2])),
                  ifelse(all(HPDinterval(mcmc(ratio[,2]))>0)|all(HPDinterval(mcmc(ratio[,2]))<0),0.05,NA))
  meccs[3,1:4] <- c(posterior.mode(mcmc(ratio[,3])),HPDinterval(mcmc(ratio[,3])), (1-length(ratiodiff[,1][ratiodiff[,1]<0])/rrows))						
  
  ## Output Verbiage -------------
  cat(paste("Here are the derived Statistics - some acos() warnings are expected out of ",rrows, "iterations","\n","\n",
            "For Ovaskainen and Angle Between Gmax, Ovaskainen's psi statistic is reported and the value in the last column is absolute", "\n",
            "Strong evidence of differences in D or Angles are evidenced by values <0.05, corresponding to 95% or more of the HPD being","\n",
            "negative; i.e. the differences between the two environments are much larger than the differences within", "\n","\n",
            "For all others, a value of 0.05 indicates signficant differences (e.g. the CI does not include 0 and all differences are positive or negative).","\n",
            "A value of NA represents non-signficance","\n","\n",
            "TCI will be 0 if there are only 3 traits","\n","\n",
            "An -> after the statistic designates a metric that compares the two environments/populations","\n","\n"))
  
  return(list(CoreStats=round(mcore,3),Vectors=mvecs,Nd_Evolv=round(Nd_evol,3),Eccentricity=round(meccs,3)))
  
}
