\name{predict.localDiscrimination}
\alias{predict.localDiscrimination}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction function for local discrimination model
}
\description{
Computed aggregated posterior probabilities for a local discrimination model.
}
\usage{
predict.LocalDiscrimination(cores, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cores}{A local discrimination object computed by localDiscrimination function.}
  \item{data}{A data.frame with the same structure as the training dataset for the localDiscrimination function.}
}

\value{
\item{weights}{Aggregated posterior probabilities from each core, weighted by the expected quality of the respective core. Each row refers to an observation, each column to a possible class.}
\item{class}{Class with the highest posterior probability from weights}
\item{weights.raw}{Unweighted, aggregated posterior probabilities of each core. Each row refers to an observation, each column to a possible class.}
\item{class.raw}{Class with the highest posterior probability from weights.raw}
}
\references{

}
\author{Thomas Ortner}

\seealso{
See Also as \code{\link{localDiscrimination}}
}


