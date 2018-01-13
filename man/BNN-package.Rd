\name{BNN-package}
\alias{BNN-package}
\docType{package}
\title{
Bayesian Neural Networks for High-Dimensional Nonlinear Variable Selection
}
\description{ 
Perform Bayesian variable selection for high-dimensional nonlinear systems. The computation can be accelerated using multiple CPUs. The package can also be used to test nonlinearity for a general regression problem. 
}
\details{
\tabular{rl}{
Package: \tab BNN\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2018-01-08\cr
License: \tab GPL-2\cr
}
The Bayesian neural network used in the package is a one-hidden layer feedforward neural network with shortcut connections.
The first module of the package is to calculate the prior probabilities assigned to the class of linear models (i.e., those 
 networks with only shortcut connections) and the class of nonlinear models; and the second module is to perform nonlinear variable selection and calculate the posterior probabilities of the classes of linear and nonlinear models. 
}
\author{
Bochao Jia, Faming Liang
Maintainer: Bochao Jia<jbc409@ufl.edu>
}
\references{
Liang, F., Li, Q. and Zhou, L. (2017) Bayesian Neural Networks for Selection of Drug Sensitive Genes.  Journal of the American Statistical Association, under revision. 

Xue, J. and Liang, F. (2017). Robust model-free feature screening for ultrahigh dimensional data.  Journal of Computational and Graphical Statistics, in press.
}
\examples{
#library(BNN)
#BNNprior(50, 1, hid_num = 3,lambda=0.025,total_iteration = 1000000, popN = 20)
}
\keyword{ package }

