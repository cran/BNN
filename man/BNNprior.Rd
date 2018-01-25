\name{BNNprior}
\alias{BNNprior}

\title{ Prior Probability of Bayesian Neural Networks
}
\description{
Calculating the prior probability of linear and nonlinear classes of BNN models.
}
\usage{
BNNprior(dimX, dimY, hid_num = 3,lambda=0.025,total_iteration=1000000,popN = 20)
}
\arguments{
  \item{ dimX }{ Dimension of the input data.}
  \item{ dimY  }{The dimension of reponse data. It is restricted to 1 in the current version of the package.}
  \item{ hid_num  }{ Number of hidden units. The default setting is 3.}
    \item{lambda}{The prior probability for each connection of the neural network being selected for the final model.
                The default setting is 0.025.}
  \item{ total_iteration  }{Number of total iterations, default of 1000,000.}
  \item{ popN  }{ Number of Markov Chains, default of 20.}
}
\value{
  \item{prob}{Prior probability assigned to the class of linear models.}
    
%% ...
}
\author{ Bochao Jia and Faming Liang}
\examples{
library(BNN)
BNNprior(50, 1, hid_num = 3, lambda=0.025, total_iteration = 10000, popN = 5)
}
\references{
Liang, F., Li, Q., and Zhou, L. (2017). Bayesian Neural Networks for Selection of Drug Sensitive Genes. Journal of the American Statistical Association.

}

\keyword{BNNprior}

