\name{Kmo.test}
\alias{Kmo.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{KMO test
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Kmo.test(corr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{corr}{corr matrix
%%     ~~Describe \code{corr} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{1001001ii
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (corr) 
{
    ifsolve <- T
    ifsolve <- tryCatch({
        r <- solve(corr)
        ifsolve <- T
    }, error = function(e) {
        massage("Matrix is not invertible.")
        return(F)
    })
    if (ifsolve == T) {
        r <- cov2cor(r)
        diag(corr) <- 0
        diag(r) <- 0
        sumr2 <- sum(r^2)
        sumcorr2 <- sum(corr^2)
        Kmo <- sumcorr2/(sumcorr2 + sumr2)
        cat("Kmo:", Kmo, sep = "")
        return(Kmo)
    }
    else return()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
