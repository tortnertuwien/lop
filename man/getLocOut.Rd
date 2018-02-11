\name{getLocOut}
\alias{getLocOut}
\title{
Get Local Outlyingness
}
\description{
Calculates the local description (getLocalDescription) for each observation of a data set. getLocalDescription is computed for each observation of data. To increase performance the package parrallel is used.
}
\usage{
getLocOut(data, k, alpha=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A data set to be evaluated for outliers. A matrix or a data frame may be handed over.}
  \item{k}{The number of observations used to define the environment by the k nearest neighbors. Only [k*alpha] of those observations are used for the definition of the projection space and for the covariance estimation used for the calculation of core distances.}
  \item{alpha}{A robustness parameter to adjust the number of oberservations used for the definition of the projection space. The optimal alpha depends on the expected number of outliers. Unless k is very small (<10), or no additional information is available, alpha should be kept at 0.5.}
}

\value{
  \item{ods}{Matrix of orthogonal distances. The columns refere to the projection, the rows to the observation.}
  \item{cds}{Matrix of core distances. The columns refere to the projection, the rows to the observation.}
  \item{ecds}{Matrix of Euclidean core distances. The columns refere to the projection, the rows to the observation.}
  \item{core}{Matrix, documenting the [k*alpha] observations used for projections (columns). 1 referes to used, 0 to not used.}
  \item{knn}{Matrix, documenting the k nearest neighbors for each observation under evaluation. 1 referes to the knn, 0 to not not part of knn.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Thomas Ortner (thomas.ortner@tuwien.ac.at)
}


\examples{
library(locout)
LocOut <- getLocOut(glass$data, k=20, alpha=0.5)
}

