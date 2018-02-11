\name{weightedLocOut}
\alias{weightedLocOut}

\title{
Weighted Outlyingness
}

\description{
Calculates the weighted oulyingness based on orthogonal weighted by inverse core distances.}

\usage{
weightedLocOut(LocOut)
}

\arguments{
\item{lo}{Object calculated by getLocOut. Must be a list, including od and cd, two matrices of orthogonal and score distances for each observation (rows) and each projection (columns).}
}

\details{}

\value{
\item{}{A univariate measure for outlyingness of each observation} }

\references{
%% ~put references to the literature/web site here ~
}
\author{
Thomas Ortner (thomas.ortner@tuwien.ac.at)
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(lop)
LocOut <- getLocOut(glass$data, k=20, alpha=0.5)
w <- weightedLocOut(LocOut)
plot(w, col=glass$group)


#Note that group one, even though considered as a consistant
#group in the original publication, holds additional structure,
#identified as outliers. For more details on the data set see
#Lemberge, P., De Raedt, I., Janssens, K. H., Wei, F., and
#Van Espen, P. J. (2000). Quantitative analysis of 16–17th
#century archaeological glass vessels using pls regression of
#epxma and μ-xrf data. Journal of Chemometrics, 14(5-6):751–763.
}
