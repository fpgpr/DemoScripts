#Pantelis has limited some of the original functions' functionality aka. THEY ARE NOT THE SAME AS BEFORE (in some cases)
#This work is only limited to cover the functions we use to produce the results in the paper
source("ICA_functions_SHOWCASE.R")

one.custom.bump <- function(u, std, len){
    #generate a gaussian bump over 100 datapoints, padded to 200
    #defaults to random positioning
    Z = rnorm(2000, mean=u, sd=std)           #generate 100k random numbers
    return(density(Z, n=len, from=0, to=1)$y) #get their histogram out
}
#================================================================================================

custom.bumps <- function(means, sds, len){
    # generate bumps from vectors of means and sds
    bumps <- sapply(1:length(means), function(x) one.custom.bump(means[x], sds[x], len)) #make the bumps
    bumps# / max(bumps) #Pantelis commented this
}
#================================================================================================

generate.custom.tree_DecemberVariant <- function(N, tree=tree,ICA.func=PCA_cubica, bumps = FALSE, theta.matrix,...){ 
    # This function is a chopped down version from the original generate.custom.tree written by David
    # It specifically uses a variant of the mixing matrix where you have noiseless correlation on the ancestral nodes
    # Outputs an object containing the original tree of size N, the bumps, and thetas used,
    # the "true/correct" mixing matrix of the sample, the mixing matrix given by ICA,
    # the distance matrix between the nodes and tips, the ICA defined basis as well as the "TRUE signals" for all the tree.

    DM <- dist.nodes(tree)      		#distance matrix between all nodes/tips in the tree
    MM <- get.mixing.matrix_NoAncNoise(theta.matrix, DM, num.bumps=ncol(bumps), whiten=FALSE)      #calls PhylogeneticNoise()
    ancestral.mixing.coefs <- MM[(N+1),]	#Set the ancestral mixing coefficients 
    signals <- bumps %*% t(MM)                  #Mix the basis signals to get all the signals in the tree

    colnames(signals) <- paste("t",colnames(DM),sep="") 

    mymeans <- rowMeans( (signals))		#Set the mean vector
    ICA <- ICA.func(signals - mymeans, ICs=nrow(theta.matrix))
						#ICA decomposion of ALL the available signals
    ICA$MM <- t(ginv(ICA$y)) %*% as.matrix(signals)
						#Mixing matrix given back by ICA

    #Remove stuff so the different reps fit in memory: 
    list(bumps=bumps, tree=tree, DM=DM, Thetas=theta.matrix,
            MM=MM, signals=signals, mean.fun = mymeans, ICA_MM = ICA$MM, ICA_sig = ICA$y)
}
#================================================================================================

get.mixing.matrix_NoAncNoise <- function(theta.matrix, DM, num.bumps, whiten = TRUE){
    # combines noise vectors to a matrix
    # theta: OU hyperparameters c(sf,l,sn) 
    # DM: cophentic distance matrix
    # noise: vector of 0's for random noise, 1's for phylogenetic noise for bumps
    # num.bumps: number of rows in the matrix
    # whiten: TRUE/FALSE z-transform Y
    sapply(1:num.bumps, function(x) make.some.noise_NotOnTheAnc(theta.matrix[x,], DM, whiten))
}
#================================================================================================

make.some.noise_NotOnTheAnc <- function(theta, DM , whiten){
    #generate noise for the tips of a tree
    #theta: OU hyperparameters c(sf,l,sn)
    #DM: cophentic distance matrix 
    #whiten: TRUE/FALSE z-transform Y
    Y = PhylogeneticNoise_NotOnTheAnc(theta, DM)  	#Generated Phylogenetic Gaussian noise for the tips 
    if(whiten) return(Whitener(Y))			#Explicitly centre to zero and make variance equal to 1
    else return(Y)
}
#================================================================================================

PhylogeneticNoise_NotOnTheAnc <- function( Theta, X){ 
#Theta: OU hyperparameters
#X : Distance matrix
# returns vector of phylogenetic related "variates"
    N = dim(X)[1]          #Number of tips  
    s_f = (  Theta[1]^2 )      #function's amplitude
    l =   (  Theta[2]   )     #characteristic lengh scale
    s_n = (  Theta[3]^2 )     #noise's amplitude
 
    N <- M <-  dim(X)[1]
    # phylogenetic GP OU VCV matrix:
    K  <- matrix( rep(0, N^2) , ncol= N)
    for ( i in 1:N){
        for ( j in i:N){ 
            K[i,j] <-  s_f * exp(-(abs( X[i,j] ))/(l)) # + s_c
        }
    }  
    K <- K + t(K)
    K <- K + diag(N) * ( .000001 - s_f) #  -s_c		#correct the double sf dose on the diagonal
    for (u in 1: ((N+1)*.5)){ K[u,u] = K[u,u]+s_n }

    VectorOfUncorrelatedRN <- rnorm(N); 		#newline
    #Q <- length( c(1: -1+(N+1)*.5 )); VectorOfUncorrelatedRN[1:Q] <- rnorm(Q);#newline
    #print(length(VectorOfUncorrelatedRN))
    return( as.numeric(t(   VectorOfUncorrelatedRN %*%   (chol(K))    )) ) #Correlated random noise
}
#================================================================================================

Whitener <- function( Coeffs){ 
#Coeffs: Vector to centre to zero and make Unit variance returns "white" vector 
    return ( (Coeffs - mean(Coeffs))/sd(Coeffs) )
}
#================================================================================================

Predictions_ForAspecificNode <- function(Y,X,Theta,new_X){  
    #Theta: OU hyperparameters
    #X : Pairwise distances of known points in the phylogeny
    #Y : Values at the tips    
    #new_X : Distance between the point of estimation and the known tips.
    # returns two-element list. 1st element : Means, 2nd element Variances
    N <- length(Y);#Number of tips   
    M <- dim(X)[1]
    #new_X : Pairwise distances of new point in the phylogeny with the respect to the existing ones. 
    
    K  <- matrix( rep(0, N^2) , ncol= N)#covariance matrix
    
    s_f = (  Theta[1]  ) #function's amplitude
    l =   (  Theta[2]  ) #characteristic lengh scale
    s_n = (  Theta[3]  ) #noise's amplitude 
        
    for ( i in 1:N){
        for ( j in i:N){ 
            K[i,j] = s_f^2 * exp(-(abs( X[i,j] ))/(l))  #+ s_c
        }
    } 
    K = K + t(K)
    K = K + diag(N) * (s_n^2 - s_f^2)# -s_c

    #Calculate the COVAR from each data point to the new one:
    K_x_xs = rep( 0,N)
    K_xs_xs = s_f^2 + s_n^2 # + s_c
    #Set Means and Vars
    Means =  0
    Vars  =  0
            
    xs_coord = new_X
    for (i in 1:N){
        K_x_xs[i] = s_f^2  * exp(  -(abs(  xs_coord[i]))/(l)  )# + s_c
    }
    Means =  K_x_xs %*% solve(K  ,Y)
    Vars  =  K_xs_xs -  K_x_xs %*% solve( K , K_x_xs)  
    return(list(Means=Means,Vars=Vars))
}
#================================================================================================

#Ancestral Node Reconstruction
Reconstructed_Curve <-function( Means= M, Vars= V, Basis = ICs, Sample_Mean = Sample_Mean ){
    #Mean Curve
    MeanCurve = t(Means) %*% Basis;
    #Standard Deviation
    Z = sqrt(rowSums(t(Basis^2) %*% diag(Vars) ))
    #3 column matrix with lower, mean and upper CI for prediction
    PredMat =  matrix(c(MeanCurve + Sample_Mean ,MeanCurve- (1.96*Z) + Sample_Mean ,  MeanCurve + (1.96*Z) + Sample_Mean ), ncol=3)
    #x11(); matplot(PredMat,type='l' ) ;
    return(PredMat)
}
#================================================================================================

Full_Curve_Estimate <- function(Y , X, thetas, new_X , Basis , Sample_Mean  = rep(0,1024)){

 P1 = Predictions_ForAspecificNode (Y= Y[1:128,1], X= X, Theta= thetas[1,], new_X= new_X)
 #print( Y[1:128,1]  )
 P2 = Predictions_ForAspecificNode (Y= Y[1:128,2], X= X, Theta= thetas[2,], new_X= new_X) 
 P3 = Predictions_ForAspecificNode (Y= Y[1:128,3], X= X, Theta= thetas[3,], new_X= new_X) 
 Ms= c( P1$Means, P2$Means, P3$Means );  Vs= c( P1$Vars, P2$Vars, P3$Vars );

 #print(Vs)
 Curve_Estimate= Reconstructed_Curve( Means= Ms, Vars= Vs, Basis =  t(Basis), Sample_Mean ) 

 return (Curve_Estimate)
}


#================================================================================================

PlotCurveEstimate <- function(Reconstructed_Curve_Obj, ylow, yhigh, extraCurve,XPLANE, strng){
    matplot(x = XPLANE,Reconstructed_Curve_Obj,ylim=c(ylow, yhigh), type= 'n', main =strng,lty =1,ylab= "y",xlab= "x")
    abline(v=(seq(0,1 ,.2)), col="lightgray", lty="dotted")
    abline(h=(seq(ylow,yhigh ,2.5)), col="lightgray", lty="dotted")
    #plot grey area polygon with samples 2STD ranges
    polygon(x=c(XPLANE[1:1024], XPLANE[1024:1]), y = c(Reconstructed_Curve_Obj[,3] ,Reconstructed_Curve_Obj[1024:1,2] ), col= 'lightgrey',border=1)
    lines(x=XPLANE, y= Reconstructed_Curve_Obj[,1], col="red", lwd = 2)
    lines(x=XPLANE, y= extraCurve, lwd = 2  )
}

#================================================================================================

PlotCurveEstimateVariationSplit <- function(Reconstructed_Curve_Obj, Phylo_Curve_Obj, ylow, yhigh, extraCurve,XPLANE, strng){
    matplot(x = XPLANE,Reconstructed_Curve_Obj,ylim=c(ylow, yhigh), type= 'n', main = strng,lty =1,ylab= "y",xlab= "x")
    abline(v=(seq(0,1 ,.2)), col="lightgray", lty="dotted")
    abline(h=(seq(ylow,yhigh ,2.5)), col="lightgray", lty="dotted")
    #plot grey area polygon with samples 2STD ranges
    #Print the whole variation
    polygon(x=c(XPLANE[1:1024], XPLANE[1024:1]), y = c(Reconstructed_Curve_Obj[,3] ,Reconstructed_Curve_Obj[1024:1,2] ), col= 'lightgrey',border=NA)
    #Print the variation only due to phylogeny
    polygon(x=c(XPLANE[1:1024], XPLANE[1024:1]), y = c(Phylo_Curve_Obj[,3] ,Phylo_Curve_Obj[1024:1,2] ), col= 'darkgrey',border=NA)
    #Put the borders on the whole variation 
    polygon(x=c(XPLANE[1:1024], XPLANE[1024:1]), y = c(Reconstructed_Curve_Obj[,3] ,Reconstructed_Curve_Obj[1024:1,2] ), col=NA,border=NA)
    lines(x=XPLANE, y= Reconstructed_Curve_Obj[,1], col="red",lwd = 2)
    lines(x=XPLANE, y= extraCurve,lwd = 2 )
}


#================================================================================================

Natural_Scale <- function( Solution ){
    #puts hyperparameter values in natural scale
    #Solution: OU hyperparamaters (in log scale) 
    # returns Solution in natural scale
    return (exp(Solution))
}

#================================================================================================

GetThetas <-function(Distances, Data){
  #This is not a generic function, it assumes that Data is 128-by-1 and Distances 128-by-128
  Ksc = Distances;
  Q1  = Data;  

  ff<-c();
  for (u in 1:100){
    SubSample <- sort(sample.int(n=128,size=100));  Ksc100 <- Ksc[SubSample,SubSample];  Q1_100 <- Q1[SubSample];
    S0 <-  uobyqa( log(runif(3)) ,fLogLik_General_Only_F,  X=Ksc100,Y=Q1_100)  
    ff<- c(ff, Natural_Scale( S0$par), S0$fval)
  }
  GG <- matrix(ff, nrow=4)
  #Sanity Check for length scales
  San <- which(GG[2,] > max(Ksc));
  GG[2,San] <- max(Ksc)
  San <- which(GG[1,]  <GG[3,]);
  GG[2,San] <- mean(GG[2,-San])

  ThetasMean= rowMeans(GG)[1:3] 
  ThetasMedian= c( median(GG[1,]), median(GG[2,]), median(GG[3,]) )  
  S<-uobyqa( log(runif(3)) ,fLogLik_General_Only_F,  X=Ksc , Y=Q1);
  ThetasStraight <- Natural_Scale(S$par); 

  return( list(ThetasMean, ThetasStraight, ThetasMedian))
}

#================================================================================================

PhylogeneticNoise <- function( Theta, X){ 
#Theta: OU hyperparameters
#X : Distance matrix
# returns vector of phylogenetic related "variates"
    N = dim(X)[1]          #Number of tips  
    s_f = (  Theta[1]^2 )      #function's amplitude
    l =   (  Theta[2]   )     #characteristic lengh scale
    s_n = (  Theta[3]^2 )     #noise's amplitude
 
    N <- M <-  dim(X)[1]
    # phylogenetic GP OU VCV matrix:
    K  <- matrix( rep(0, N^2) , ncol= N)
    for ( i in 1:N){
        for ( j in i:N){ 
            K[i,j] <-  s_f * exp(-(abs( X[i,j] ))/(l)) # + s_c
        }
    }  
    K <- K + t(K)
    K <- K + diag(N) * (s_n - s_f) #  -s_c
    
    return( as.numeric(t(  rnorm(N)  %*%   (chol(K))    )) ) #Correlated random noise
}
