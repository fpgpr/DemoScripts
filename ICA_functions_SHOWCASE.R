cubica34 <- function(x){
#
# This is just a port of the original cubica34.m in R all credits to the original authors
# mistakes are most probably from Pantelis Hadjipantelis (C:) . 2011-04-17
#
##################################################################################
#
# CubICA (IMPROVED CUMULANT BASED ICA-ALGORITHM)
#
# This algorithm performes ICA by diagonalization of third- and
# fourth-order cumulants simultaneously.
#
#  [R,y]=cubica34(x)
#
# - x is and NxP matrix of observations 
#     (N: Number of components; P: Number of datapoints(samplepoints)) 
# - R is an NxN matrix such that u=R*x, and u has 
#   (approximately) independent components.
# - y is an NxP matrix of independent components
#  
# This algorithm does exactly (1+round(sqrt(N)) sweeps.
# 
# Ref: T. Blaschke and L. Wiskott, "An Improved Cumulant Based
# Method for Independent Component Analysis", Proc. ICANN-2002,
# Madrid, Spain, Aug. 27-30.
#
# questions, remarks, improvements, problems to: t.blaschke@biologie.hu-berlin.de.
#
# Copyright : Tobias Blaschke, t.blaschke@biologie.hu-berlin.de.
#
# 2002-02-22
#
#
# Last change:2003-05-19 
#  
################################################################################## 

N = dim(x)[1] 
P = dim(x)[2] 
Q=  diag(N);
resolution=0.001;

# centering and whitening

M = rowMeans(x)
x = as.matrix(x)- M;   
E3 = eigen(x %*% t(x) /P ); D = E3$values; V=-E3$vectors;
V <- V[,ncol(V):1]; D = as.real(rev(D));
W = diag((D)^(-.5)) %*% t(V)  
y=W%*%x; 
# start rotating

for (t in 1:(1+round(sqrt(N)))){ 
    for (i in 1:(N-1) ){
        for (j in (i+1):N) {
            
            #calculating the new cumulants
            Rs =  c(i,j);
            u=y[Rs,]; 
                
            sq=u^2; 
            
            sq1=t(sq[1,])
            sq2=t(sq[2,])
            u1=(u[1,])
            u2=(u[2,]) 
            
            C111=(sq1%*%u1)/P;
            C112=(sq1%*%u2)/P;
            C122=(sq2%*%u1)/P;
            C222=(sq2%*%u2)/P;
            
            C1111=(sq1%*%t(sq1))/P -3;
            C1112=(sq1*t(u1))%*%u2 /P;
            C1122=sq1%*%t(sq2)/P -1;
            C1222=(sq2*t(u2))%*%u1 /P;
            C2222=sq2%*%t(sq2)/P -3;
            
            # coefficients
            
            c_34=(1/6)*(1/8)*(3*(C111^2+C222^2)-9*(C112^2+C122^2)-6*(C111*C122+C112*C222));
            
            c_44=(1/24)*(1/16)*(7*(C1111^2+C2222^2)-16*(C1112^2+C1222^2)-12*(C1111*C1122+C1122*C2222)-36*C1122^2-32*C1112*C1222-2*C1111*C2222);
            
            s_34=(1/6)*(1/4)*(6*(C111*C112-C122*C222));
            
            s_44=(1/24)*(1/32)*(56*(C1111*C1112-C1222*C2222)+48*(C1112*C1122-C1122*C1222)+8*(C1111*C1222-C1112*C2222));
            
            c_48=(1/24)*(1/64)*(1*(C1111^2+C2222^2)-16*(C1112^2+C1222^2)-12*(C1111*C1122+C1122*C2222)+36*C1122^2+32*C1112*C1222+2*C1111*C2222);
            
            s_48=(1/24)*(1/64)*(8*(C1111*C1112-C1222*C2222)-48*(C1112*C1122-C1122*C1222)-8*(C1111*C1222-C1112*C2222));
            
            phi_4=-atan2(s_34+s_44,c_34+c_44);
            phi_8=-atan2(s_48,c_48);
            
            B_4=sqrt((c_34+c_44)^2+(s_34+s_44)^2);
            B_8=sqrt(c_48^2+s_48^2);
            
            #calculating the angle
            
            approx=-phi_4/4-(pi/2)*trunc(-phi_4/pi);
            
            intervall=seq( (approx-pi/8),(approx+pi/8) ,resolution )
            
            psi_34=B_8*cos(8*intervall+phi_8)+B_4*cos(4*intervall+phi_4);
            
            value= max(psi_34);
            index= which.max(psi_34);
            
            phi_max=intervall[index]; 
            
            #Givens-rotation-matrix Q_ij
            
            Q_ij=diag(N);

            c=cos(phi_max);
            s=sin(phi_max);
            
            Q_ij[i,j]=s;
            Q_ij[j,i]=-s;
            Q_ij[i,i]=c;
            Q_ij[j,j]=c;
            
            Q=Q_ij%*%Q;

            # rotating y
            
            y[c(i, j),]= matrix( c(c,-s,s,c), nrow=2)%*%u;
        
        } #j
    } #i
} #t
    
R=Q%*%W;

return (list(R = R, y = y, Q=Q, W=W, x=x, u=u))
}


reorder.by.max <- function(M){
    #returns M with columns rearranged by absolute max function values
    M[,order(sapply(1:ncol(M), function(x) which(abs(M[,x])==max(abs(M[,x])))))]
}

PPCA_cubica <- function(signals, nPcs=4){
    #output of PPCA + cubica. Columns of y reordered by abs(max(col))
    print("running pPCA + CUBICA")
    prob.pca <- pca(signals, method="ppca", nPcs=nPcs, scale="uv")
    prob.ica <- cubica34(t(prob.pca@scores))
    prob.ica$y <- t(reorder.by.max(t(prob.ica$y)))
    #mm <- ica$y %*% as.matrix(signals)
    #mm
    prob.ica
}

PCA_cubica <- function(signals, ICs, ...){
    #output of PCA + cubica. Columns of y reordered by abs(max(col))
    print("running eigen PCA + CUBICA")
    PCA <- eigen(cov(t(signals)))
    ica <- cubica34(t(PCA$vectors[,1:ICs]))
    ica$y <- t(reorder.by.max(t(ica$y)))
    ica$Evals <- PCA$values
    ica$Evectors <- PCA$vectors
    ica$y <- flip.ICs(ica$y, ...)
    ica
}

PPCA_fastica <- function(signals, nPcs=4){
    #output of PPCA + cubica. Columns of y reordered by abs(max(col))
    print("running pPCA + CUBICA")
    ppca <- pca(signals, method="ppca", nPcs=nPcs, scale="uv", ...)
    ica <- cubica34(t(ppca@scores))
    ica$y <- t(reorder.by.max(t(ica$y)))
    #mm <- ica$y %*% as.matrix(signals)
    #mm
    ica
}

PPCA_fastica <- function(signals, ...){
    print("running pPCA + fastICA")
    ppca <- pca(signals, method="ppca", nPcs=4, scale="uv", ...)
    ica <- fastICA(ppca@scores, n.comp = 4, row.norm = TRUE)
    ica$y <- t(reorder.by.max(ica$S))
    
}

flip.vector <- function(x, flip = TRUE, zero = TRUE){
    #flips a vector if the maximum value is negative
    # aligns lhs to zero
    if(flip) if(max(x) != max(abs(x))) x <- -x
    if(zero) x <- x + (0-x[1])
    x
}

flip.ICs <- function(ICs, flip.vectors = TRUE, zero.intercept = TRUE){
    #apply flip.vector to all columns in an IC matrix 
    ICs <- t(ICs)
    flipped <-sapply(1:ncol(ICs), 
        function(x) flip.vector(ICs[,x],flip=flip.vectors, zero=zero.intercept))
    t(flipped)
}


PPCA_cubica.flip <- function(signals, nPcs){
    # PPCA + cubica, flipped and baseline brought to zero
    ica <- PPCA_cubica(signals=signals, nPcs=nPcs)
    ica$y <- flip.ICs(ica$y)
    ica
}

library(caTools);
Normalize_To_Area1 <- function(y,x){
	#crude normalization to area 1 using the trapezoid rule.
	#Assuming N by L matrix 
	L = dim(y)[2] ; N = dim(y)[1] ;  
	Normalized_y= y;
	for (i in 1:N ){Normalized_y[i,] =  y[i,] /  trapz(y= y[i,], x=x) }
    return ( Normalized_y  )
}

