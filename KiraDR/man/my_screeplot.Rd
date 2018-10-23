\name{my_screeplot}
\alias{my_screeplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Screeplot
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
my_screeplot(pca, n, title = "Screeplot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pca}{a pca matrix
%%     ~~Describe \code{pca} here~~
}
  \item{n}{number of factors
%%     ~~Describe \code{n} here~~
}
  \item{title}{the title of picture
%%     ~~Describe \code{title} here~~
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
function (pca, n, title = "Screeplot") 
{
    data_var <- apply(pca[, 1:n], 2, var)
    data_plot <- data.frame(comp = c(1:n), var = data_var)
    ggplot(data_plot, aes(x = comp, y = var, group = 1)) + geom_line(size = 1) + 
        geom_point(size = 3) + theme(axis.title.x = element_text(size = 25, 
        face = "bold.italic", vjust = 0.5, hjust = 0.5)) + theme(axis.title.y = element_text(size = 25, 
        face = "bold.italic", vjust = 0.5, hjust = 0.5)) + theme(axis.text.x = element_text(size = 15, 
        vjust = 0.5, hjust = 0.5)) + theme(axis.text.y = element_text(size = 15, 
        vjust = 0.5, hjust = 0.5)) + scale_x_continuous(breaks = seq(0, 
        n, n\%/\%5)) + ggtitle(title) + theme(plot.title = element_text(size = 30))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
