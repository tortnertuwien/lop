\name{guidedProjections}
\alias{guidedProjections}

\title{
Guided Projections Data Transformation
}

\description{
Computes the orthogonal and score distances and OSD for guided projections.}

\usage{
guidedProjections(data, k=10, osd=c("OD", "SD", "OD*SD")
}

\arguments{
\item{data}{A data matrix of n observations and p variables.}
\item{k=10}{The number of observations used to define the projections.}
\item{OSD}{Used as similarity measure for the replacement of observations and for the comparison of observations. OSD can be set to "OD" for orthogonal distances, "SD" for score distances and "OD*SD" for the product of orthogonal and score distances.}
}

\details{Computes the orthogonal and score distances and OSD for a series of projections, providing a description of the datastructure. Each projection is corresponds to the projection onto a space spanned by k observations. The observations get sequencially replaced based on a measure of similarity. Guided projections are capeable of revealing group structure and outliers in a dataset. The transformation is especially helpful for the following situations:
* inhomogenous groups of observations (e.g. different covariance structure)
* groups located in different subspaces}

\value{
A list-object of class "guidedprojections"" containing the following matrices:
\item{OD}{A matrix of orthogonal distances for each observation (columns) to each projection space (rows)}
\item{SD}{A matrix of score distances for each observation (columns) in each projection space (rows)}
\item{OSD}{A matrix measures of similarity for each observation (columns) and projection (rows)}
}

\references{
Guided projections for analysing the structure of high-dimensional data
https://arxiv.org/abs/1702.06790
}
\author{
Thomas Ortner (thomas.ortner@tuwien.ac.at)
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(lop)
#might take a minute
gp <- guidedProjections(glass$data)
plot(gp, col=glass$group)

#For details on the data set see
#Lemberge, P., De Raedt, I., Janssens, K. H., Wei, F., and
#Van Espen, P. J. (2000). Quantitative analysis of 16–17th
#century archaeological glass vessels using pls regression of
#epxma and μ-xrf data. Journal of Chemometrics, 14(5-6):751–763.
}
