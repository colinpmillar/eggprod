\name{get_survey_area}
\alias{get_survey_area}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get survey area
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_survey_area(year = c("1998", "1998_extended", "1995"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{year}{
%%     ~~Describe \code{year} here~~
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
\author{
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
function (year = c("1998", "1998_extended", "1995")) 
{
    year <- match.arg(year)
    if (year == "1995") {
        list(nmin = 17L, nmax = 44L, longmin = c(1, 1, 1, 1, 
            1.5, 2, 2.5, 4, 4.5, 5, 6, 6.5, 7, 7.5, 8, 10, 10, 
            10, 10, 10, 10, 9, 8, 8, 8, 8, 8, 8), longmax = c(3, 
            5, 5, 6, 6.5, 9.5, 11, 12, 13, 14, 14, 14, 14, 14, 
            16, 16, 16, 16, 16, 16, 16, 13.5, 11.5, 11, 10.5, 
            10.5, 10.5, 10.5))
    }
    else if (year == "1998") {
        list(nmin = 17L, nmax = 45L, longmin = c(1, 1, 1, 1, 
            1.5, 2, 2.5, 4, 4.5, 5, 6, 6.5, 7, 7.5, 8, 10, 10, 
            10, 10, 10, 10, 9, 8, 8, 8, 8, 8, 8, 8), longmax = c(3, 
            5, 5, 6.5, 10, 10, 14, 14, 14, 14, 14, 14, 14, 16, 
            16, 16, 16, 16, 16, 16, 16, 13.5, 11.5, 11, 10.5, 
            10.5, 11, 11, 11))
    }
    else if (year == "1998_extended") {
        list(nmin = 17L, nmax = 46L, longmin = c(1, 1, 1, 1, 
            1.5, 2, 2.5, 4, 4.5, 5, 6, 6.5, 7, 7.5, 8, 10, 10, 
            10, 10, 10, 10, 9, 8, 8, 8, 8, 8, 8, 8, 8), longmax = c(3, 
            5, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
            20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
            20, 20, 20))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
