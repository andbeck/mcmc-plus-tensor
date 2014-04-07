########################################################################################
# psb2() [plotsubspace2] function by Beckerman and Robinson 2011
#
# The function plots either G matrices; E matrices; P matrices (i.e. population plasticity) 
# In all cases, if No. Traits = 3, it is the empirical space.  
# If traits > 3 it is the first three principal components of the matrix
# When n > 3, axes are labeled by eigenvectors (i.e. loadings) and 
#
# The function requires arguments of 
# a) two models
# b) the method - G or P
# c) the number of traits
# d) graphical details about shading and axes labels.
#
# G plots centred covariance hulls for the G-structure from a MCMCglmm model
# Could be modified to plot E-matrix : centred Environmental covariance hulls for the R-structure from a MCMCglmm model

# --Both of these have error checking to ensure that the model has random effects.
# --This is basically Hadfield (MCMCglmm package) plotsubspace() function with some extra outputs
# --and designed to take models rather than matrices.

# P plots the Residual covariance structure from a MCMCglmm model with no random effect
# --the Hulls in the P method are displaced by the distance between the traits in multivariate space
# --A small sphere marks the centroid of each hull and a line connects them.
# --The P method describes the change in MEAN and VARIANCE of the traits.
# --This is a modified plotsubspace() function to allow the offsetting and add centroids
# 
# Two additional functions are provided that are internal to our psb2:
# --ellisoid3d
# --collect.matrices
########################################################################################

###############################################################
################### ELLIPSOID3d() #############################
############# INTERNAL FUNCTION in plotsubspace() #############
############# Actually from the rgl method ####################
###############################################################

ellipsoid3d<-function(rx=1,ry=1,rz=1,n=30,ctr=c(0,0,0), qmesh=FALSE,trans = par3d("userMatrix"),...) {
	
	if(require(rgl)==FALSE){stop("rgl not loaded")}

#-----------------------------------------------------------------#
## This function is taken from Daniel Adlers and Duncan Murdochs ##
##                         rgl package                           ##
##                     as used by Hadfield                    	 ##
#-----------------------------------------------------------------#

if (missing(trans) && !rgl.cur()) trans <- diag(4)

degvec <- seq(0,2*pi,length=n)

ecoord2 <- function(p) {
	c(rx*cos(p[1])*sin(p[2]),ry*sin(p[1])*sin(p[2]),rz*cos(p[2])) }
	
	v <- apply(expand.grid(degvec,degvec),1,ecoord2)
	if (qmesh) v <- rbind(v,rep(1,ncol(v))) ## homogeneous
	
	e <- expand.grid(1:(n-1),1:n)
	i1 <- apply(e,1,function(z)z[1]+n*(z[2]-1))
	i2 <- i1+1
	i3 <- (i1+n-1) %% n^2 + 1
	i4 <- (i2+n-1) %% n^2 + 1
	i <- rbind(i1,i2,i4,i3)
		if (!qmesh)
			quads3d(v[1,i],v[2,i],v[3,i],...)
		else return(rotate3d(qmesh3d(v,i,material=...),matrix=trans))
}


###############################################################################################
###############################################################################################
# collect.matrices: Function to collect relevant matrices from the MCMCglmm models
###############################################################################################
###############################################################################################

collect.matrices<-function(model,no.traits,shrink=FALSE){
	
	if (!inherits(model, "MCMCglmm")) 
        stop("use only with \"MCMCglmm\" objects")
	
	# get dimension of Gcov
	GCtest<-length(as.numeric(summary(model)$Gcovariances[,1]))
	
	#Gmatrix calculation	
	ifelse(GCtest!=0, GCmat<-matrix(posterior.mode(model$VCV[,1:GCtest]),sqrt(GCtest),sqrt(GCtest),byrow=TRUE), GCmat <- 0)
	
	# Residual P (E) matrix calculation
	ifelse(GCtest==0, RCmat<-matrix(posterior.mode(model$VCV),no.traits,no.traits,byrow=TRUE), RCmat<-0) # we ignore this if method = G.... only valid when method = P
	
	# Phen values for centroids
	P<-as.numeric(summary(model)$solutions[,1])

	# make generic trait names
	ll<-length(as.data.frame(model$Sol))
	nn<-paste("trait",1:no.traits)
	
	return(list(GenCov=GCmat,RCov=RCmat,Locate=P,Names=nn))
	}


########################################################################################
## PLOTSUBSPACE 2 - psb2
########################################################################################

psb2 <- function (model1,model2,no.traits,method=c("G","P"),corr = FALSE, shadeCA = TRUE, shadeCB = TRUE, 
    axes.lab = FALSE,...) 
{
	if (require(rgl) == FALSE)
        stop("rgl not loaded")
	if (!inherits(model1, "MCMCglmm")) 
        stop("use only with \"MCMCglmm\" objects")
	
	# use collect matrices function to get all cov mats and location paramters
	# specify shrinkage use here

	e1mat<-collect.matrices(model1,no.traits,shrink=FALSE)
	e2mat<-collect.matrices(model2,no.traits,shrink=FALSE)		

	# Follow Jarrod methods for subspace plotting - G or R or P or I   
	#-------------------------------------------------------------------------#
	## This function is modified from Jarrod Hadfields plotsubspace function ##
	##                         MCMCglmm package                              ##
	#-------------------------------------------------------------------------#
	
	#----------------------------- G Method method ---------------------------#
	if(method=="G"){
		if(is.null(dim(e1mat$GenCov)))
    		stop("Method G - No Gmatrix calculated in model")
    		
	if(method=="G"){
		CA<-e1mat$GenCov
				dimnames(CA)<-list(e1mat$Names,e1mat$Names)
		CB<-e2mat$GenCov
			dimnames(CB)<-list(e1mat$Names,e1mat$Names)
		}
	
	if(method=="E"){
		CA<-e1mat$RCov
			dimnames(CA)<-list(e1mat$Names,e1mat$Names)
		CB<-e1mat$RCov
			dimnames(CB)<-list(e1mat$Names,e1mat$Names)
		}
		
	Avec <- eigen(CA)$vectors[, 1:3]
	Aval <- eigen(CA)$values
	Aval[Aval<0]<-0										# fudge in the event there is a negative eigenvalue to ALLOW plotting; BE CAREFUL.
	propA <- round(100 * sum(Aval[1:3])/sum(Aval), 0)
   	Aval <- 2 * sqrt(Aval[1:3])
	
	clear3d()
	s1 <- ellipsoid3d(Aval[1], Aval[2], Aval[3], qmesh = TRUE, 
        trans = diag(4))
    if (dim(CA)[1] == 3) {
        s1 <- rotate3d(s1, matrix = t(Avec))
    }
    if (shadeCA == TRUE) {
        shade3d(s1, col = "red")
    }
    else {
        wire3d(s1, col = "red")
    }
    if (is.null(CB) == FALSE) {
        if (dim(CA)[1] == 3) {
            Bvec <- eigen(CB)$vectors
            Bval <- 2 * sqrt(eigen(CB)$values)
            propB<-round(100 * sum(Bval)/sum(eigen(CB)$values),0)
        }
        else {
            Bproj <- t(Avec) %*% CB %*% Avec
            Bvec <- eigen(Bproj)$vectors
            Bval <- eigen(Bproj)$values
            Bval[Bval<0]<-0										# fudge in the event there is a negative eigenvalue to ALLOW plotting; BE CAREFUL.
            propB <- round(100 * sum(Bval)/sum(eigen(CB)$values), 0)
            Bval <- 2 * sqrt(Bval)
        }
        s2 <- ellipsoid3d(Bval[1], Bval[2], Bval[3], qmesh = TRUE, 
            trans = diag(4))
        s2 <- rotate3d(s2, matrix = t(Bvec))
        if (shadeCB == TRUE) {
            shade3d(s2, col = "blue")
        }
        else {
            wire3d(s2, col = "blue")
        }
    }
   if (dim(CA)[1] < 4) {
        xvec <- e1mat[1]
        yvec <- e1mat[2]
        zvec <- e1mat[3]
    		}
   else {
        xvec <- paste("[", paste(round(Avec[, 1], 2), collapse = ","), 
            "]", sep = "")
        yvec <- paste("[", paste(round(Avec[, 2], 2), collapse = ","), 
            "]", sep = "")
        zvec <- paste("[", paste(round(Avec[, 3], 2), collapse = ","), 
            "]", sep = "")
    		}
    nameA <- as.character(substitute(CA))
    nameA <- nameA[length(nameA)]
    if (is.null(CB)) {
        nameB <- NULL
    }
    else {
        nameB <- as.character(substitute(CB))
        nameB <- nameB[length(nameB)]
    }
    axes3d()
    
    # Title and Axes Labels
    if (axes.lab == TRUE) {
        if(method=="G"){title3d(main="Genetic Covariance Matrices",
        	xlab = paste("PC1: ",xvec), ylab = paste("PC2: ",yvec), zlab = paste("PC3: ", zvec))}
        if(method=="E"){title3d(main="Genetic Covariance Matrices",
        	paste("PC1: ",xvec), ylab = paste("PC2: ",yvec), zlab = paste("PC3: ",zvec))}
        }   
    
	cat("Loadings","\n")
	load<-matrix(round(Avec,2),length(e1mat$Names),3,dimnames=list(c(e1mat$Names),c("PC1","PC2","PC3")))
	print(load)
	cat("\n")
	cat("Proportion of Variation in PC1-3","\n")
	props<-matrix(c(propA,propB),1,2,dimnames=list(c("Proportion"),c("Env 1","Env2")))
	print(props)	
   }
#----------------------- P Method ----------------------------------------------------------------------------------------#
	if(method=="P"){
	
	if(!is.null(dim(e1mat$GenCov)))
		stop("Do not use P method with models estimating Genetic Covariance")
	
		CA<-e1mat$RCov
			dimnames(CA)<-list(e1mat$Names,e1mat$Names)
		CB<-e2mat$RCov
			dimnames(CB)<-list(e2mat$Names,e2mat$names)
	
	Avec <- eigen(CA)$vectors[, 1:3]
	Aval <- eigen(CA)$values
	Aval[Aval<0]<-0										# fudge in the event there is a negative eigenvalue to ALLOW plotting; BE CAREFUL.
	propA <- round(100 * sum(Aval[1:3])/sum(Aval), 0)
	Aval <- 2 * sqrt(Aval[1:3])
	
	clear3d()
	
	s1 <- ellipsoid3d(Aval[1], Aval[2], Aval[3], qmesh = TRUE,trans = diag(4))
	 if (dim(CA)[1] == 3){
        s1 <- rotate3d(s1, matrix = t(Avec))}
    		# ADDED TRANSLATION HERE TO MOVE SUBSPACE - require use of wire to see inside
    		wire3d(translate3d(s1,e1mat$Locate[1],e1mat$Locate[2],e1mat$Locate[3]),col = "red")
        
        if (is.null(CB) == FALSE) {
        
        if (dim(CA)[1] == 3) {
            Bvec <- eigen(CB)$vectors
            Bval <- 2 * sqrt(eigen(CB)$values)
            propB<-round(100 * sum(Bval)/sum(eigen(CB)$values),0)
        }
        else {
            Bproj <- t(Avec) %*% CB %*% Avec
            Bvec <- eigen(Bproj)$vectors
            Bval <- eigen(Bproj)$values
            Bval[Bval<0]<-0										# fudge in the event there is a negative eigenvalue to ALLOW plotting; BE CAREFUL.
            propB <- round(100 * sum(Bval)/sum(eigen(CB)$values), 0)
            Bval <- 2 * sqrt(Bval)
        }
        
        s2 <- ellipsoid3d(Bval[1], Bval[2], Bval[3], qmesh = TRUE, trans = diag(4))
        s2 <- rotate3d(s2, matrix = t(Bvec))
    		# ADDED TRANSLATION HERE TO MOVE SUBSPACE - require wire to see inside
    		wire3d(translate3d(s2,e2mat$Locate[1],e2mat$Locate[2],e2mat$Locate[3]),col = "blue")
    		
  		}
  		
  	if (dim(CA)[1] < 4) {
        xvec <- e1mat[1]
        yvec <- e1mat[2]
        zvec <- e1mat[3]
    		}
    else {
        xvec <- paste("[", paste(round(Avec[, 1], 2), collapse = ","), 
            "]", sep = "")
        yvec <- paste("[", paste(round(Avec[, 2], 2), collapse = ","), 
            "]", sep = "")
        zvec <- paste("[", paste(round(Avec[, 3], 2), collapse = ","), 
            "]", sep = "")
    		}
    nameA <- as.character(substitute(CA))
    nameA <- nameA[length(nameA)]
    
    if (is.null(CB)) {
        nameB <- NULL
    }
    else {
        nameB <- as.character(substitute(CB))
        nameB <- nameB[length(nameB)]
    }
    axes3d(expand=1)
 
    # calulate HPDinterval on difference of posterior location estimates for each trait
    
    ptit<-	round(HPDinterval(model1$Sol-model2$Sol),2)
    
    # Title and Axes Labels
    if (axes.lab == TRUE) {
        title3d(main="Phenotypic Covariance Matrices",
        	paste("PC1: ",xvec), ylab = paste("PC2: ",yvec), zlab = paste("PC3: ",zvec))
    }
	# add centroids and line here to mark plasticity?
  	rgl.spheres(x=c(e1mat$Locate[1],e2mat$Locate[1]),
  				y=c(e1mat$Locate[2],e2mat$Locate[2]),
  				z=c(e1mat$Locate[3],e2mat$Locate[3]),
  				color=c("red","blue"),radius=0.25,alpha=0.5)
	
	rgl.lines(x=c(e1mat$Locate[1],e2mat$Locate[1]),
				y=c(e1mat$Locate[2],e2mat$Locate[2]),	
				z=c(e1mat$Locate[3],e2mat$Locate[3]),
				color="black",lwd=3)
	
	# Return Plasticity Estimates
	if(length(e1mat$Locate)>3) cat("Warning - Centroids May be wrong in plot","\n\n")
	cat("HPDinterval of Differences in Trait Locations (Do they move?)","\n")
	cat("--- All Pos or All Neg by Row = Significant Plasticity","\n")
	print(ptit)
	cat("\n")
	cat("Loadings","\n")
	load<-matrix(round(Avec,2),length(e1mat$Names),3,dimnames=list(c(e1mat$Names),c("PC1","PC2","PC3")))
	print(load)
	cat("\n")
	cat("Proportion of Variation in PC1-3","\n")
	props<-matrix(c(propA,propB),1,2,dimnames=list(c("Proportion"),c("Env 1","Env2")))
	print(props)				
	}
}

########################################################################################
#----------------------- END PLOT METHODS ---------------------------------------------#
########################################################################################

