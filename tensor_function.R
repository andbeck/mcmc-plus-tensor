#### -----------------------------------------------------------------------------------------
#### Begin tensor.stats FUNCTION
#### By Matt Robinson and Andrew Beckerman
####
## Now Hosted on BitBucket (7.4.14)
## NOTE CURRENTLY EXPECT MODELS FIT TO RETURN JOINT POSTERIOR WITH LENGTH 1000
## ENSURE that your models have nitts (!) and burnin's set to deliver you 1000
#### -----------------------------------------------------------------------------------------

tensor.stats<-function(models,no.traits,no.trt){
	if (require(MCMCglmm) == FALSE)
		stop("MCMCglmm not loaded")
	if (require(mvtnorm) == FALSE)
        stop("mvtnorm not loaded")
	if (!inherits(models[[1]], "MCMCglmm")) 
        stop("use only with \"MCMCglmm\" objects")
    if(no.trt!=length(models))
        stop("Missing a model for some treatments")
	if(dim(models[[1]]$VCV)[1]!=1000)
		stop("please re-run your models to get a joint posterior of 1000 samples")
	
	n<-no.traits
	m<-no.trt
	d<-(n*(n+1)/2)+1 # tensor dimension (when is it +1????)
		
	wrk.set<-vector("list", m)
	
	for(i in 1:m){
		wrk.set[[i]]<-data.frame(models[[i]]$VCV[,1:(n)^2])
	}

# Collection Zones
	sum_eigvals<-matrix(0,1000,d)
	sum_projmat<-vector("list",1000)
	sum_eigtenvecs<-vector("list",1000)
	sum_eigtenvals<-vector("list",1000)

# Set up tensor shite lists
	proj_mat<-vector("list",(d-1)) # list has n*(n+1)/2) matrices (rows of projmat) each m = trt columns
	for (i in 1:(d-1)){
		proj_mat[[i]]<-matrix(0,1000,m)
	}

	eigten_vec<-vector("list",d) # list has n*(n+1)/2)+1 matrices each n = trait columns 
	for (i in 1:(d-1)){
		eigten_vec[[i]]<-matrix(0,1000,n)
	}

	eigten_vals<-vector("list",d) # list has n*(n+1)/2)+1 matrices each n = trait columns 
	for (i in 1:(d-1)){
		eigten_vals[[i]]<-matrix(0,1000,n)
	}

# MASTER LOOPING OVER 1000
	for (z in 1:1000){	

########################################################################################################
## 1. Stack the matrices from each MCMC model into one file where ncol = number of traits and nrow = n traits * i environments
		Gs<-list()
			for(i in 1:m){
			Gs[[i]]<-matrix(as.numeric(wrk.set[[i]][z,][,1:(n)^2]),n,n)
			}
	
		Gmatrices<-matrix(unlist(Gs),m*n,n,byrow=TRUE)

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
			count<-1
		  	for (i in 1:(n-1)){
    			for (j in (i+1):n){
		      		covcov[k,count+n]=sqrt(2)*cov(Gmatrices[seq(k,n*m,n),k],Gmatrices[seq(i,n*m,n),j])
      				covcov[count+n,k]=sqrt(2)*cov(Gmatrices[seq(k,n*m,n),k],Gmatrices[seq(i,n*m,n),j])
				count=count+1
				}
			}
		}

# lower right quadrant of S
		countx<-1
		county<-1
		for (k in 1:(n-1)){
			for (l in (k+1):n){
				for (i in 1:(n-1)){
					for (j in (i+1):n){
						covcov[countx+n,county+n]=2*cov(Gmatrices[seq(k,n*m,n),l],Gmatrices[seq(i,n*m,n),j])
						covcov[county+n,countx+n]=2*cov(Gmatrices[seq(k,n*m,n),l],Gmatrices[seq(i,n*m,n),j])
						county=county+1}
        			}
			countx=countx+1
			county=1
			}
		}

		covcov<-as.data.frame(covcov)

########################################################################################################
## 3. Take eigenvectors and eigenvalues of covariance tensor S, 
## and then for each eigenvector of S construct and eigentensor,
## then take eigenvectors and eigenvalues of each eigentensor

		eigvecs<-eigen(covcov)$vectors
		eigvals<-eigen(covcov)$values
		eigten<-matrix(nrow=n*n*(n+1)/2,ncol=n)

		for (eval in 1:(n*(n+1)/2)){ 								#number of eigenvalues of covariance tensor
			for (i in 1:n){ 										#constructing an n*n eigentensor from each eigenvector of S matrix
    			eigten[i+n*(eval-1),i]<-eigvecs[i,eval] 			#this loop calculates the diagonal
			}

			count<-1
				for (i in 1:(n-1)){
					for (j in (i+1):n){  						eigten[i+n*(eval-1),j]<-1/sqrt(2)*eigvecs[count+n,eval]
						eigten[j+n*(eval-1),i]<-1/sqrt(2)*eigvecs[count+n,eval]
					count<-count+1
					}
				}
		}

# storage matrices of eigentensor eigenvalues and eigenvectors
		eigtenvecs<-matrix(nrow=n*n*(n+1)/2, ncol=n)
		eigtenvals<-matrix(nrow=n*n*(n+1)/2, ncol=1)

		for (i in 1:(n*(n+1)/2)){
			eigtenvecs[((i-1)*n+1):(i*n),]<-t(eigen(eigten[((i-1)*n+1):(i*n),])$vectors) # eigenvectors in rows as in Excel
			eigtenvals[((i-1)*n+1):(i*n),1]<-eigen(eigten[((i-1)*n+1):(i*n),])$values
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

#  ADJUST 1-2 level collections here and at start to be lists..... to make generic

		# Summary of eigenvalues
		sum_eigvals[z,1:(d-1)]<-eigvals
		sum_eigvals[z,d]<-sum(eigvals)

		# Projection Matrix Collection
		sum_projmat[[z]]<-projmat

		# summary vectors and values from proj matrix
		sum_eigtenvecs[[z]]<-eigtenvecs
		sum_eigtenvals[[z]]<-eigtenvals

		# allocate details
		for(i in 1:(d-1)){
			proj_mat[[i]][z,]<-projmat[i,]
			}

		for(i in 1:(d-1)){
			eigten_vec[[i]][z,]<-eigtenvecs[i,]
			}
		
		code<-rep(1:(n*m), each=n)
		for(i in 1:(d-1)){
			eigten_vals[[i]][z,]<-eigtenvals[code==i]
			}
	}

	tens.out<-list(
		sum_eigvals,
		sum_projmat,
		sum_eigtenvecs,
		sum_eigtenvals,
		proj_mat,
		eigten_vec,
		eigten_vals,
		tr.names=rownames(summary(models[[1]])$solutions))
	
	return(tens.out)
}
#### -----------------------------------------------------------------------------------------
#### END tensor.stats FUNCTION
#### -----------------------------------------------------------------------------------------


#### -----------------------------------------------------------------------------------------
#### Begin tensor.analyse FUNCTION
#### By Matt Robinson and Andrew Beckerman
#### -----------------------------------------------------------------------------------------

tensor.analyse<-function(dd,how.many.tensors=1, no.traits = 2){
	n<-no.traits
	if(length(dd)!=8)
	stop("tensor stats function failed to produce correct info")
cat("\n***********************************\n")
cat("*** Core Tensor Analysis Output ***\n")
cat("***********************************\n")	
cat("\n --- 95% CI of Variance Explained by each tensor \n")
	tmp<-c()
	for(i in 1:how.many.tensors){
		cat("\n Axis ",i," \n")
		tmp[i]<-posterior.mode(as.mcmc(dd[[1]][,i]/dd[[1]][,((n*(n+1)/2)+1)]))
		print(
			cbind(
				MODE = posterior.mode(as.mcmc(dd[[1]][,i]/dd[[1]][,((n*(n+1)/2)+1)])),
				CI = HPDinterval(as.mcmc(dd[[1]][,i]/dd[[1]][,((n*(n+1)/2)+1)]))
				)
			)
	}
cat("\nThe number of significant tensors is ~ where the index of this cumsum = 1 \n")
cat("Consider adding another tensor dimension if it has not reached 1 \n")

print(cumsum(tmp))

cat("\n *For the next two derived stats, know your traits:\n")
print(dd[[8]])

cat("\n --- For the significant eigentensors, the Directional Change in Variance of each trait for each vector")
	for(i in 1:how.many.tensors){
		cat("\n Axis ",i," \n")
		print(
			cbind(
			MODE = posterior.mode(as.mcmc(dd[[7]][[i]][,1:no.traits])), 
			CI = HPDinterval(as.mcmc(dd[[7]][[i]][,1:no.traits]))
				)
			)
	}	
	
cat("\n\n --- The Proportion of Variance of each eigentensor explianed by each trait \n")
for(i in 1:how.many.tensors){
		cat("\n Axis ",i," \n")
		tmp<-dd[[7]][[i]]
		tdim<-dim(tmp)[2]
		pvar<- matrix(0,1000,no.traits)
		# calc tot var across all traits
		# should be 1 for 1st and <1 for next
		tot<-apply(tmp,1,function(x) sum(x^2))
		# for each trait
			for(j in 1:1000){
				for(k in 1:no.traits){
					pvar[j,k]<-tmp[j,k]^2/tot[j]	
					}
				}
		print(
			cbind(
				MODE = posterior.mode(as.mcmc(pvar[,1:no.traits])), 
				CI = HPDinterval(as.mcmc(pvar[,1:no.traits]))
				)
			)
	}

cat("\n *For the next stats, recall the order of your treat/envs in the model list* \n")

cat("\n --- Genetic variance of each of GEI for each env")
	for(i in 1:how.many.tensors){
		cat("\n Axis ",i," \n")
		tdim<-dim(dd[[5]][[i]])[2]
		print(
			cbind(
				Mode = posterior.mode(as.mcmc(sqrt(dd[[5]][[i]][,1:tdim]^2))),
				CI = HPDinterval(as.mcmc(sqrt(dd[[5]][[i]][,1:tdim]^2)))
				)
			)
	}

cat("\n\n *** Extra Output ***")
cat("\n 95% CI of Eigenvalues of the Eigentensors \n")
	print(posterior.mode(as.mcmc(dd[[1]][,1:((n*(n+1)/2)+1)])))
	print(HPDinterval(as.mcmc(dd[[1]][,1:((n*(n+1)/2)+1)])))

}