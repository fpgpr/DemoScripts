#==================================
#Getting the tree to use
#==================================
rm(list=ls())
library(ape) 	#to use various phylo-related goodies
library(car)	#to use Box-Cox transform
library(geiger) #to use fitContinuous() to get λ
library(moments)#to use kurtosis and skewness
library(minqa)	#to use uobyqa for the quadratic approximation in optimization 
	
N = 128; 	#Number of extant taxa

rr=22;
set.seed(98); 
RTree <- rtree(N, br= exp( rnorm(N)*sqrt(2.22) - 2.45)  );  
		#Produce the tree we could like  (set seed to fix rtree()'s branching events)
#pdf("Figure_Tree.pdf", width=5, height=5)
x11(); 
plot.phylo(RTree, show.tip.label=F)
#nodelabels(129:255,cex =.333)
#dev.off();
		
source("Function_Collection_SHOWCASE.R"); #
	

 
#1 also gives nice ICA vs PCA results.

set.seed(rr);

#==================================
#Getting the basis signals to use
#==================================
#Specify signal parameters (Length, means, "width of bell")
sig.len <- 1024				
sig.means <- c(200, 500, 700 )/sig.len
sig.sds <- c(60, 120, 70 )/sig.len
#Produce signals.
sig.matrix <- custom.bumps(sig.means, sig.sds, sig.len)

tree <- RTree; # (max(cophenetic(tree)) == 8.870363)

#Set the thetas
thetas <-matrix(c(2.5, 0.75* max(cophenetic(tree)),.5, 		 	 
                  0,   max(cophenetic(tree)),  1,   	 
                  1.5, 0.25*max(cophenetic(tree)),.5),	 
               ncol=3, byrow=TRUE) # matrix of thetas for all bumps

#Set the seed here again, you might want to sample some nice decompositions.

#==================================
#Getting the tips signals to use
#==================================
#Generate the signals throughout the tree so we can use them for validation afterwards
xbase <- generate.custom.tree_DecemberVariant (N= 128, tree=tree, ICA.func=PCA_cubica, bumps = sig.matrix, theta.matrix=thetas);
#xbase <- generate.custom.tree_whole(N= 128, tree=tree, ICA.func=PCA_cubica, bumps = sig.matrix, theta.matrix=thetas) 
XPLANE = seq(0,1,length.out=1024);

#Generate the plot with the signals on their own so you can copy paste to the bigger file. 
if(1==2){
pdf("Nodes_For_PlotR.pdf", width=10, height=50);
par(mfrow=c(4,1)) 
plot( y=xbase$signals[,67],  x= XPLANE, ylim=c(-20,20), lwd =14.5,type='l', ann=FALSE, xaxt ='n', yaxt='n', bty='n');
plot( y=xbase$signals[,5], x= XPLANE, ylim=c(-20,20), lwd = 14.5,type='l', ann=FALSE, xaxt ='n', yaxt='n', bty='n');
plot( y=xbase$signals[,105], x= XPLANE, ylim=c(-20,20), lwd = 14.5,type='l', ann=FALSE, xaxt ='n', yaxt='n', bty='n');
plot( y=xbase$signals[,33],  x= XPLANE, ylim=c(-20,20), lwd = 14.5,type='l', ann=FALSE, xaxt ='n', yaxt='n', bty='n');
dev.off();

pdf("Nodes_For_PlotR_internal.pdf", width=10, height=25);
par(mfrow=c(3,1)) 
plot( y=xbase$signals[,233],  x= XPLANE, ylim=c(-20,20), lwd = 14.25,type='l', ann=FALSE, xaxt ='n', yaxt='n', bty='n');
plot( y=xbase$signals[,129], x= XPLANE, ylim=c(-20,20), lwd = 14.25,type='l', ann=FALSE, xaxt ='n', yaxt='n', bty='n');
plot( y=xbase$signals[,151],  x= XPLANE, ylim=c(-20,20), lwd = 14.25,type='l', ann=FALSE, xaxt ='n', yaxt='n', bty='n');
dev.off();
}


#==================================
#Getting the 4 suplots figure.
#==================================
#pdf("Figure_IPCA.pdf")
x11(); 
source("GettingTheIPCAFigure_SHOWCASE.R") #This might be slightly different than the one the paper
#dev.off();

 
 
if(1==1){

#plot(mlbench.smiley(500, .1,0.05), xlim=c(-1.65,1.65))


#==================================
#Getting the Normality Plots of the Mixing Matrices
#==================================

#source("GettingTheNormalityPlots_SHOWCASE.R")

#==================================
#Getting the 3 subplots estimates fig.
#==================================

#Compile the MLE fuctions.
source("cxxFunctions_OU_Generic_SHOWCASE.cxx")
wert =2;
set.seed(wert)

Ksc <- cophenetic(RTree)
Nthetas <-thetas;
Nthetas[1,] <- GetThetas(Data=xbase$ICA_MM[1,1:N], Ksc)[[1]]
if (Nthetas[1,1] < Nthetas[1,3]) { Nthetas[1,] = c(0.00001,1, sd( xbase$ICA_MM[1,1:N])) }
Nthetas[2,] <- GetThetas(Data=xbase$ICA_MM[2,1:N], Ksc)[[1]]
if (Nthetas[2,1] < Nthetas[2,3]) { Nthetas[2,] = c(0.00001,1, sd( xbase$ICA_MM[2,1:N])) }
Nthetas[3,] <- GetThetas(Data=xbase$ICA_MM[3,1:N], Ksc)[[1]]
if (Nthetas[3,1] < Nthetas[3,3]) { Nthetas[3,] = c(0.00001,1, sd( xbase$ICA_MM[3,1:N])) }

#Get Estimated Thetas
 Estimated_thetas = Nthetas; #These are the stuff that get feed in the estimates figure plotter!
print((Nthetas-thetas)/thetas)
#This script plots pdf stuff "inside" it.
x11(); 
source("GettingTheEstimatesFigure_SHOWCASE.R") 
 

}


