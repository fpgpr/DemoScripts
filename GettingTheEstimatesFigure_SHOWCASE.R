
E_thetas_OLD =  Estimated_thetas 
E_thetas =  Estimated_thetas; E_thetas[,3]= 0.00001;
E_thetasClean <- Estimated_thetas; E_thetasClean[,1:3] <- 0.00001;
E_thetasClean_OLD<- Estimated_thetas; E_thetasClean_OLD[,1:2] <- 0.00001;
if(1==3){

pdf("Fig4_131_13Dec_EstimatedThetas.pdf",width=14, height=5)
RootEstimate = Full_Curve_Estimate(Y= t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetas, new_X = dist.nodes(tree)[129,1:128], Basis =t(xbase$ICA_sig) )
NodeEstimate=Full_Curve_Estimate(Y=t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetas, new_X = dist.nodes(tree)[131,1:128], Basis = t(xbase$ICA_sig))
TipEstimate =Full_Curve_Estimate(Y=t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetas_OLD, new_X = dist.nodes(tree)[033,1:128], Basis = t(xbase$ICA_sig) ) 

RootEstimateC = Full_Curve_Estimate(Y= t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetasClean, new_X = dist.nodes(tree)[129,1:128], Basis =t(xbase$ICA_sig) )
NodeEstimateC=Full_Curve_Estimate(Y=t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetasClean, new_X = dist.nodes(tree)[131,1:128], Basis = t(xbase$ICA_sig))
TipEstimateC =Full_Curve_Estimate(Y=t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetasClean_OLD, new_X = dist.nodes(tree)[033,1:128], Basis = t(xbase$ICA_sig) )  

RootEstimateC =RootEstimate -RootEstimateC
NodeEstimateC =NodeEstimate -NodeEstimateC
TipEstimateC  = TipEstimate -TipEstimateC
YUPPER = max( abs(c(RootEstimate, TipEstimate, NodeEstimate)))
print(YUPPER)

par(mfrow=c(1,3)); 
PlotCurveEstimateVariationSplit(XPLANE=XPLANE, ylow=-YUPPER,yhigh=+YUPPER,extraCurve=xbase$signals[,129],Reconstructed_Curve_Obj=RootEstimate,Phylo_Curve_Obj =RootEstimateC, strng ="Root Estimate" )
PlotCurveEstimateVariationSplit(XPLANE=XPLANE, ylow=-YUPPER,yhigh=+YUPPER,extraCurve=xbase$signals[,131],Reconstructed_Curve_Obj= NodeEstimate,Phylo_Curve_Obj = NodeEstimateC, strng ="Node Estimate")
PlotCurveEstimateVariationSplit(XPLANE=XPLANE, ylow=-YUPPER,yhigh=+YUPPER,extraCurve=xbase$signals[,033], Reconstructed_Curve_Obj= TipEstimate, Phylo_Curve_Obj = TipEstimateC, strng ="Tip Estimate")
dev.off()



pdf("Fig4_233_13Dec_EstimatedThetas.pdf",width=14, height=5)
RootEstimate = Full_Curve_Estimate(Y= t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetas, new_X = dist.nodes(tree)[129,1:128], Basis =t(xbase$ICA_sig) )
NodeEstimate=Full_Curve_Estimate(Y=t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetas, new_X = dist.nodes(tree)[233,1:128], Basis = t(xbase$ICA_sig))
TipEstimate =Full_Curve_Estimate(Y=t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetas_OLD, new_X = dist.nodes(tree)[033,1:128], Basis = t(xbase$ICA_sig) ) 

RootEstimateC = Full_Curve_Estimate(Y= t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetasClean, new_X = dist.nodes(tree)[129,1:128], Basis =t(xbase$ICA_sig) )
NodeEstimateC=Full_Curve_Estimate(Y=t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetasClean, new_X = dist.nodes(tree)[233,1:128], Basis = t(xbase$ICA_sig))
TipEstimateC =Full_Curve_Estimate(Y=t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetasClean_OLD, new_X = dist.nodes(tree)[033,1:128], Basis = t(xbase$ICA_sig) ) 

RootEstimateC =RootEstimate -RootEstimateC
NodeEstimateC =NodeEstimate -NodeEstimateC
TipEstimateC  = TipEstimate -TipEstimateC
YUPPER = max( abs(c(RootEstimate, TipEstimate, NodeEstimate)))
print(YUPPER)
par(mfrow=c(1,3)); 
PlotCurveEstimateVariationSplit(XPLANE=XPLANE, ylow=-YUPPER,yhigh=+YUPPER,extraCurve=xbase$signals[,129],Reconstructed_Curve_Obj=RootEstimate,Phylo_Curve_Obj =RootEstimateC, strng ="Root Estimate" )
PlotCurveEstimateVariationSplit(XPLANE=XPLANE, ylow=-YUPPER,yhigh=+YUPPER,extraCurve=xbase$signals[,233],Reconstructed_Curve_Obj= NodeEstimate,Phylo_Curve_Obj = NodeEstimateC, strng ="Node Estimate")
PlotCurveEstimateVariationSplit(XPLANE=XPLANE, ylow=-YUPPER,yhigh=+YUPPER,extraCurve=xbase$signals[,033], Reconstructed_Curve_Obj= TipEstimate, Phylo_Curve_Obj = TipEstimateC, strng ="Tip Estimate")
dev.off()
}



#pdf("Fig4_151_13Dec_EstimatedThetas.pdf",width=14, height=5)
RootEstimate = Full_Curve_Estimate(Y= t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetas, new_X = dist.nodes(tree)[129,1:128], Basis =t(xbase$ICA_sig) )
NodeEstimate=Full_Curve_Estimate(Y=t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetas, new_X = dist.nodes(tree)[151,1:128], Basis = t(xbase$ICA_sig))
TipEstimate =Full_Curve_Estimate(Y=t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetas_OLD, new_X = dist.nodes(tree)[033,1:128], Basis = t(xbase$ICA_sig) ) 

RootEstimateC = Full_Curve_Estimate(Y= t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetasClean, new_X = dist.nodes(tree)[129,1:128], Basis =t(xbase$ICA_sig) )
NodeEstimateC=Full_Curve_Estimate(Y=t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetasClean, new_X = dist.nodes(tree)[151,1:128], Basis = t(xbase$ICA_sig))
TipEstimateC =Full_Curve_Estimate(Y=t(xbase$ICA_MM), X=cophenetic(tree), thetas=E_thetasClean_OLD, new_X = dist.nodes(tree)[033,1:128], Basis = t(xbase$ICA_sig) ) 

RootEstimateC =RootEstimate -RootEstimateC
NodeEstimateC =NodeEstimate -NodeEstimateC
TipEstimateC  = TipEstimate -TipEstimateC
YUPPER = max( abs(c(RootEstimate, TipEstimate, NodeEstimate)))
print(YUPPER)
par(mfrow=c(1,3)); 
PlotCurveEstimateVariationSplit(XPLANE=XPLANE, ylow=-YUPPER,yhigh=+YUPPER,extraCurve=xbase$signals[,129],Reconstructed_Curve_Obj=RootEstimate,Phylo_Curve_Obj =RootEstimateC, strng ="Root Estimate" )
PlotCurveEstimateVariationSplit(XPLANE=XPLANE, ylow=-YUPPER,yhigh=+YUPPER,extraCurve=xbase$signals[,151],Reconstructed_Curve_Obj= NodeEstimate,Phylo_Curve_Obj = NodeEstimateC, strng ="Node Estimate")
PlotCurveEstimateVariationSplit(XPLANE=XPLANE, ylow=-YUPPER,yhigh=+YUPPER,extraCurve=xbase$signals[,033], Reconstructed_Curve_Obj= TipEstimate, Phylo_Curve_Obj = TipEstimateC, strng ="Tip Estimate")
#dev.off()


