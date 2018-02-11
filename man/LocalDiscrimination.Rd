\name{localDiscrimination}
\alias{localDiscrimination}

\title{
Local Discrimination based discrimination model
}

\description{
Uses class-specific cores to compute local proections. The core distances and the orthogonal distance of each local projection are used to compute an LDA model. The models are then aggregated to an overall model.
}

\usage{
localDiscrimination(data, class, k, dist=NULL)
}

\arguments{
\item{data}{A data matrix containing the training data for the local projection based discrimination model. The observations are stored row-wise.}
\item{class}{The class of observations, stored as a factor variable.}
\item{k}{A tung parameter determining the degree of locality to be used for the local projections. It can be set to any natural number between 3 and the smallest group size-1. A good starting point is 50\% of the smallest group size.}
}

\details{}
A vector of lists of all considered local projections:ld[[1]]
\value{
\item{core}{The indices of observations used for the core.}
\item{group}{The class of the observation initializing the local projection.}
\item{projection}{A projection matrix used to project obersvations onto the core space. Computed using an SVD.}
\item{vf}{Variable Filter. Variables with no variation within the core observations are removed as they cannot be re-scaled.}
\item{location}{Estimation of location.}
\item{scale}{Estimation of scatter.}
\item{core.location}{Core based estimation of location.}
\item{core.scale}{Core based estimation of scatter.}
\item{frequency}{The frequency of the core occuring.}
\item{quality}{A class specific description of the quality of the LDA model.}
\item{cv.posterior}{Matrix of oosterior probabilities for each observation for each class computey by cross validation.}
\item{lda}{The local LDA model of class "lda"}
\item{df.core}{The centered and scaled core representation of the training observations.}
}

\references{
%% ~put references to the literature/web site here ~
Ortner, T., Hoffmann, I., Filzmoser, P., Zaharieva, M., Breiteneder, C., and Brodinova, S., ”Multigroup discrimination based on weighted local projections. Submitted to Journal of Computational and Graphical Statistics.
}
\author{
Thomas Ortner (thomas.ortner@tuwien.ac.at)
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
See Also \code{\link{predict.localDiscrimination}, \link{plot.localDiscrimination}}
}
\examples{
library(lop)

x <- rrcovHD::olitos
data <- x[,1:25]
cl <- as.factor( x[,26] )
ld <- localDiscrimination(data, cl, k=7)
p0 <- predict(ld, data)
table(p0$class, cl)

lpplot(p0$weights, cl, combination=c(3,2))
lpplot(p0$weights, cl)

}
