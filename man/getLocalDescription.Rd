\name{getLocalDescription}
\alias{getLocalDescription}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get the local description of an observation.}
\description{
Uses the k nearest neighbours of observation x. From those observations the most dense [k*alpha] observations are used to provide a projection onto the [k*alpha] oberservations. The orthogonal distance to this projection space and the score distance within this projection space, based on the covariance of the [k*alpha] observations are calculated for each observation from the dataset.
}
\usage{
getLocalDescription(data, center, dist=NULL, k=20, alpha=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A dataset, with n rows and p variables. Any matrix or data frame can be handed over.
}

  \item{center}{
An index of data, which defines the observation initiating the projection process.
}
  \item{dist}{
A distance matrix of data. If no distance matrix is provided, it will be calculated using dist(data). The parameter is used in order to make the calculation more efficient.
}
\item{k=20}{
the number of nearest neighbours considered. k further influences the dimension of the projection space [k*alpha] at the same time.
}
\item{alpha}{
The robustness parameter for the estimation of the covariance matrix. It further influences the dimension of the projection space.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{od}{A vector of orthogonal distances of each observation with respect to the space spanned by the core observations.}
\item{cd}{A vector of core distances of each observation with respect to the covariance structure defined by the core observations.}
\item{ecd}{A vector of Euclidean core distances of each observation.}
\item{core}{The indices of [k*alpha] observations used to span the projection space and estimate the covariance matrix.}
\item{knn}{The k nearest neighbours of the center observation.}
\item{center}{The observation initiating the local description.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Thomas Ortner (thomas.ortner@tuwien.ac.at)
}

\examples{
library(lop)
locDesc <- getLocalDescription(glass$data, center=1, k=20, alpha=0.5)
pch <- rep(1, nrow(glass$data))
pch[locDesc$core] <- 4
plot(locDesc$od ~ locDesc$cd, col=glass$group, pch=pch)
}
