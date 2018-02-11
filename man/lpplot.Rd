\name{lpplot}
\alias{lpplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting posterior probabilities for multigroup discrimination results}
\description{
}
\usage{
lpplot(probabilities, class, combination=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probabilities}{A matrix of posterior probabilities from a classification result. Not nessesarily from local projection based discrimination.}
  \item{class}{A factor variable classifying the observations. Either the ground truth or the classification result.}
  \item{combination}{A vector of length two or not provided at all. If provided, the plot will be reduced to the visualization of the separation of the indicated groups. If not provided, all possible combinations will be visualized in form of a matrix.}
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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

