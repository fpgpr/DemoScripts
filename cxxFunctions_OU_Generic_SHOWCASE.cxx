library(RcppEigen)
library(inline)

fLogLik_General_inR <-function(Thetas, Y, X){
#This is only done for the sake of checking the correctness of other functions.
N = length(Y);
s_f2 = exp(2*Thetas[1])
l = exp(Thetas[2])
s_n2 = exp(2*Thetas[3])
rr=0;

K = matrix( rep(0,N*N), nrow=N)
K_x_x= K;
for(i in 1:N){
  for(j in i:N){
    rr= X[i,j]
    K[i,j] = s_f2 * exp(- rr/ l)
  }
}
K_x_x = t(K)+ K ;
K = K_x_x + diag(N)*(s_n2 - s_f2 +.000001);
A1 =  -.5* t(Y) %*% solve(K,Y);
A2 =  sum(log(diag( chol(K))) )  

return( -( A1 - A2 - (N/2)* log(2*pi) ))
}


fLogLik_General <- cxxfunction(signature(ThetaI = "vector",YY = "vector", XX= "matrix"),
    '
#include <cmath>  	// to use sqrt and exp
using Eigen::Map;   	// to map input variable to an existing array of data
using Eigen::MatrixXd;  // to use MatrixXd
using Eigen::VectorXd;  // to use VectorXd
using Eigen::ArrayXd;   // to use ArrayXd
using Eigen::LLT;   	// to do the LLT decomposition
using Eigen::Lower; 	// to get the lower triangular view
const Map<VectorXd> Theta(Rcpp::as<Map<VectorXd> > (ThetaI));   //Map vector ThetaI to matrixXd Theta
const Map<VectorXd> Y(Rcpp::as<Map<VectorXd> > (YY));       	//Map vector YY to vectorXd Y
const Map<MatrixXd> X(Rcpp::as<Map<MatrixXd> > (XX));       	//Map matrix QQ to matrixXd X
//const Map<VectorXd> B(Rcpp::as<Map<VectorXd> > (BB));       	//Map vector BB to vectorXd B

using namespace std;
  
    int N= Y.size();  //number of points
    double s_f2, l, s_n2;     //hyperparameters / function amplitude, characteristic length, noise amplitude, s_c 
    ArrayXd dF = ArrayXd::Constant( Theta.size(),0);

    s_f2 = exp( 2.0* Theta(0) ); //exponentiate to make sure they are positive
    l   = exp(      Theta(1) );
    s_n2 = exp( 2.0* Theta(2) ); 
    double rr = .0;

    MatrixXd K = MatrixXd::Zero(N,N) ;        //Covariance K  
    MatrixXd K_x_x = MatrixXd::Zero(N,N) ;    //Covariance K helper
    for (int i=0; i<N; i++){
      for (int j=i; j<N; j++){
        rr = (X(i,j)); 
        K(i,j) = s_f2 * exp(-( rr )/l) ;
      }
    } 

    K_x_x = K.transpose() + K;                    // Built full matrix
    K = K_x_x + MatrixXd::Identity(N,N) * ( s_n2 - s_f2 +.000001);     // Add Noise //took out  -s_c
    LLT<MatrixXd> llt_K(K); //compute the Cholesky decomposition of the matrix
    double logLikelihood = 0.;

    if (llt_K.info()==Eigen::Success) {       // if the Cholesky decomposition exists.
      VectorXd alpha = llt_K.solve(Y);                   //Get alpha = inv(K)Y
      double A1 = 0.5*(Y.transpose()*alpha).value();         //Fit term
     double A2 = llt_K.matrixLLT().diagonal().array().log().sum();  //Complexity term
    logLikelihood = ( A1 + A2 + (N*.5)*log(2*M_PI));       //Final log-likelihood

     }
    else { 
      return (List::create(Named("F") =1234567890. ));
    }
   return (List::create(Named("F") = logLikelihood));',
plugin = "RcppEigen") 



Log_Likelihood_General <- function(Theta, Y, X){ 
# Theta : Hyperparamaters 
# Y : Trait values at tips
# X : Distance matrix
# returns logLikelihood
    storage.mode(X) <- "double";        #Save stuff as doubles
    storage.mode(Y) <- "double";
    storage.mode(Theta) <- "double"; 
 
    return ( fLogLik_General(Theta, Y, X));   #LogLikelihood
}  

fLogLik_General_Only_F  <- function(Theta, Y, X){ 
# Theta : Hyperparamaters 
# Y : Trait values at tips
# X : Distance matrix
# returns logLikelihood 
    storage.mode(X) <- "double";        #Save stuff as doubles
    storage.mode(Y) <- "double";
    storage.mode(Theta) <- "double";
 
    return (fLogLik_General(Theta, Y, X)$F);     #LogLikelihood value
} 

