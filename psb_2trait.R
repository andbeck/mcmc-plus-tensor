## Now Hosted on GitHub

###############################################################################################
###############################################################################################
# collect.matrices: Function to collect relevant matrices from the MCMCglmm models
###############################################################################################
###############################################################################################

collect.matrices<-function(model,no.traits=2){
	
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

psb_2trait <- function (model1,model2,no.traits=2,method=c("G","P"),
	corr = FALSE, shadeCA = TRUE, shadeCB = TRUE, 
    axes.lab = FALSE, xlims=NULL, ylims=NULL, xlab="Trait1", ylab="Trait2",...) 
{
		if (require(ellipse) == FALSE)
        stop("ellipse library not loaded")
		
		if (!inherits(model1, "MCMCglmm")) 
        stop("use only with \"MCMCglmm\" objects")
		
	
	# use collect matrices function to get all cov mats and location paramters
	# specify shrinkage use here

	e1mat<-collect.matrices(model1)
	e2mat<-collect.matrices(model2)		

	# Follow Jarrod methods for subspace plotting - G or R or P or I   
	#-------------------------------------------------------------------------#
	## This function is modified from Jarrod Hadfields plotsubspace function ##
	##                         MCMCglmm package                              ##
	#-------------------------------------------------------------------------#
	
	#----------------------------- G Method method ---------------------------#
	if(method=="G"){
		if(is.null(dim(e1mat$GenCov)))
    		stop("Method G inappropriate - No Gmatrix calculated in model")
    	
    	if (length(as.numeric(summary(model1)$Gcovariances[,1]))<2) 
        	stop("Random effects not speficied to estiamte VCV")
    		
    		
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
	
	es1<-eigen(CA)
	es2<-eigen(CB)
	varA<-sum(es1$values)
	varB<-sum(es2$values)
	
	plot(ellipse(CA),type="l",asp=1)
		
		Val <- es1$values
		Vec <- es1$vectors 
		v1 <- sqrt(Val[1]*qchisq(0.95,2))*Vec[,1]
		v2 <- sqrt(Val[2]*qchisq(0.95,2))*Vec[,2]
		segments(-v1[1],-v1[2],v1[1],v1[2]) 
		segments(-v2[1],-v2[2],v2[1],v2[2]) 	
		
	points(ellipse(CB),type="l",col="grey40",lty=2)
		Val <- es2$values
		Vec <- es2$vectors 
		v1 <- sqrt(Val[1]*qchisq(0.95,2))*Vec[,1]
		v2 <- sqrt(Val[2]*qchisq(0.95,2))*Vec[,2]
		segments(-v1[1],-v1[2],v1[1],v1[2],col="grey40",lty=2) 
		segments(-v2[1],-v2[2],v2[1],v2[2],col="grey40",lty=2) 	
	
	cat("\n")
	cat("Variance in Each Environment (sum Eigenvalues)","\n")
	props<-matrix(c(varA,varB),1,2,dimnames=list(c("Proportion"),c("Env 1","Env2")))
	print(props)	
   }#----------------------- P Method ----------------------------------------------------------------------------------------#
	if(method=="P"){
	
	if(!is.null(dim(e1mat$GenCov)))
		stop("Do not use P method with models estimating Genetic Covariance")
	if (require(ellipse) == FALSE)
        stop("ellipse library not loaded")
	if (!inherits(model1, "MCMCglmm")) 
        stop("use only with \"MCMCglmm\" objects")
	if (length(as.numeric(summary(model1)$Rcovariances[,1]))<2) 
        stop("RCov not speficied to estiamte VCV")
		
	CA<-e1mat$RCov
		dimnames(CA)<-list(e1mat$Names,e1mat$Names)
	CB<-e2mat$RCov
		dimnames(CB)<-list(e2mat$Names,e2mat$names)
	
	ctrA<-e1mat$Locate
	ctrB<-e2mat$Locate
	
	es1<-eigen(CA)
	es2<-eigen(CB)
	varA<-sum(es1$values)
	varB<-sum(es2$values)
	
	
	plot(ellipse(CA,centre=ctrA),type="l",asp=1, xlim=xlims, ylim=ylims, xlab=xlab, ylab=ylab)
	points(ctrA[1],ctrA[2],pch=19)
		
		Val <- es1$values
		Vec <- es1$vectors 
		v1 <- (sqrt(Val[1]*qchisq(0.95,2))*Vec[,1])
		v2 <- (sqrt(Val[2]*qchisq(0.95,2))*Vec[,2])
		segments(-v1[1]+ctrA[1],-v1[2]+ctrA[2],v1[1]+ctrA[1],v1[2]+ctrA[2]) 
		segments(-v2[1]+ctrA[1],-v2[2]+ctrA[2],v2[1]+ctrA[1],v2[2]+ctrA[2]) 	
		
	points(ellipse(CB,centre=ctrB),type="l",col="grey40",lty=2)
	points(ctrB[1],ctrB[2],pch=19,col="grey40")
	
		Val <- es2$values
		Vec <- es2$vectors 
		v1 <- (sqrt(Val[1]*qchisq(0.95,2))*Vec[,1])
		v2 <- (sqrt(Val[2]*qchisq(0.95,2))*Vec[,2])
		segments(-v1[1]+ctrB[1],-v1[2]+ctrB[2],v1[1]+ctrB[1],v1[2]+ctrB[2],col="grey40") 
		segments(-v2[1]+ctrB[1],-v2[2]+ctrB[2],v2[1]+ctrB[1],v2[2]+ctrB[2],col="grey40") 
	
	
	}
}

########################################################################################
#----------------------- END PLOT METHODS ---------------------------------------------#
########################################################################################

