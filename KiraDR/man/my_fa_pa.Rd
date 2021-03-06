\name{my_fa_pa}
\alias{my_fa_pa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Principal Factor Method
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
my_fa_pa(corr, nfactor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{corr}{corr matrix
%%     ~~Describe \code{corr} here~~
}
  \item{nfactor}{number of factors
%%     ~~Describe \code{nfactor} here~~
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
function (corr, nfactor) 
{
    p <- nrow(corr)
    diag_corr <- diag(corr)
    sum_rank <- sum(diag_corr)
    rowname <- colnames(corr)
    colname <- paste("PA", 1:nfactor, sep = "")
    loading <- matrix(0, nrow = p, ncol = nfactor, dimnames = list(rowname, 
        colname))
    k <- 1
    kmax <- 50
    if (qr(corr)$rank == p) 
        sigma <- diag(1/solve(corr))
    else sigma <- 1 - apply(abs(corr - diag(1, nrow = p)), 2, 
        max)
    h <- diag_corr - sigma
    repeat {
        diag(corr) <- h
        h1 <- h
        eig <- eigen(corr)
        for (i in 1:nfactor) loading[, i] <- sqrt(eig$values[i]) * 
            eig$vectors[, i]
        h <- diag(loading \%*\% t(loading))
        if ((sqrt(sum((h - h1)^2)) < 1e-04) | k == kmax) 
            (break)()
        k <- k + 1
    }
    rowname <- c("SS loadings", "Proportion Var", "Cumulative Var")
    var_sum <- matrix(0, nrow = 3, ncol = nfactor, dimnames = list(rowname, 
        colname))
    for (i in 1:nfactor) {
        var_sum[1, i] <- sum(loading[, i]^2)
        var_sum[2, i] <- var_sum[1, i]/sum_rank
        var_sum[3, i] <- sum(var_sum[1, 1:i])/sum_rank
    }
    method <- c("Principal Factor Method")
    list(method = method, loadings = loading, var_sum = var_sum, 
        sp = diag_corr - h, iterative = k)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
