
par(mfrow=c(2,2)) 

#==============
#Subplot 1
#==============
#print the original bumps
matplot(x = XPLANE, xbase$bumps, main = "Original Basis", ylab= "y",  lwd = 1.5, type='l', lty = 1:3, col='black', xlab= "x", ylim=c(-4,7)) 
#Make xy-y axis
abline(v=(seq(0,1 ,.2)), col="lightgray", lty="dotted");abline(h=(seq( -4,7 ,1)), col="lightgray", lty="dotted")
 
#==============
#Subplot 2
#==============
#Set the signals we know:
SignalsAtTheTips = xbase$signals[,1:128]
#Calculate the standard statistics to plot
StDevCurve  = apply(SignalsAtTheTips,1, sd); 
MeanCurve   =  rowMeans(SignalsAtTheTips)
FourCurves  = (xbase$signals[,c(67,5,105,33)]) #same curves as the one shown on the tree
#Prepare dimensions but do not print plot 
matplot(x = XPLANE ,FourCurves,  ylim=c( -2.05* max(StDevCurve ) + min( MeanCurve) , 2.05 * max(StDevCurve + max(MeanCurve)) ), type= 'n', main = "Mixed sample", lty =1, ylab= "y", xlab= "x")
abline(v=(seq(0,1 ,.2)), col="lightgray", lty="dotted")
abline(h=(seq(-40,40 ,5)), col="lightgray", lty="dotted")
#plot grey area polygon with samples 2STD ranges
polygon(x=c(XPLANE, XPLANE)  , y = c(-2.00*StDevCurve + rowMeans(SignalsAtTheTips) , +2.00*StDevCurve+ rowMeans(SignalsAtTheTips)), col= 'lightgrey',border=1)
lines(x=XPLANE, y= rowMeans(SignalsAtTheTips), col='red')
#add sample lines
lines(x = XPLANE , y= FourCurves[,1],lwd = 1.5, type='l',lty =1)
lines(x = XPLANE , y= FourCurves[,2],lwd = 1.5, type='l',lty =2)
lines(x = XPLANE , y= FourCurves[,3],lwd = 1.5, type='l',lty =3)
lines(x = XPLANE , y= FourCurves[,4],lwd = 1.5, type='l',lty =4) 

#==============
#Subploat 3
#==============
#Calculate the sample PCs 
SampleEigenvectors <- eigen(cov(t(SignalsAtTheTips)))$vectors[,1:3]
matplot(x =XPLANE ,60*  t(flip.ICs(t(SampleEigenvectors))) , main = "Principal components", ylab= "y",    lwd = 1.5, type='l', lty = 1:3, col='black', xlab= "x", ylim=c(-4,7)) #Pantelis changed this line from 1:4 to 1:3
abline(v=(seq(0,1 ,.2)), col="lightgray", lty="dotted")
abline(h=(seq(-4,7 ,1)), col="lightgray", lty="dotted")
SampleIndepComponents <- cubica34(t(SampleEigenvectors))

#==============
#Subploat 4
#==============
#Calculate the sample ICs 
ICs_To_Plot = Normalize_To_Area1 (flip.ICs(SampleIndepComponents$y), x= XPLANE)
matplot(x =  XPLANE , t( ICs_To_Plot[c(3,2,1),] ), main = "Independent\nPrincipal components" , ylab= "y",    lwd = 1.5, type='l', lty = 1:3, col='black', xlab= "x", ylim=c(-4,7))
abline(v=(seq(0,1 ,.2)), col="lightgray", lty="dotted")
abline(h=(seq(-4,7 ,1)), col="lightgray", lty="dotted")   
