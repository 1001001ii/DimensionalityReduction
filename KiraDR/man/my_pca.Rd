\name{my_pca}
\alias{my_pca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PCA
%%  ~~function to do ... ~~
}
\description{PCA
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
my_pca(data, cor = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data frame
%%     ~~Describe \code{data} here~~
}
  \item{cor}{use cor or not
%%     ~~Describe \code{cor} here~~
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
function (data, cor = F) 
{
    data <- scale(data, center = T, scale = F)
    if (cor == F) {
        sigma <- cov(data)
        eigenvector.sigma <- eigen(sigma)$vectors
        pca.sigma <- as.matrix(data) \%*\% eigenvector.sigma
        return(pca.sigma)
    }
    if (cor == T) {
        error <- F
        error <- tryCatch(r <- cor(data), warning = function(e) {
            message("Warning:The standard deviation is zero, and the correlation coefficient matrix can not be calculated, so use the covariance matrix to solve the principal component.")
            return(T)
        })
        if (error == F) {
            eigenvector.r <- eigen(r)$vectors
            pca.r <- as.matrix(data) \%*\% eigenvector.r
            return(pca.r)
        }
        else {
            sigma <- cov(data)
            eigenvector.sigma <- eigen(sigma)$vectors
            pca.sigma <- as.matrix(data) \%*\% eigenvector.sigma
            return(pca.sigma)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
