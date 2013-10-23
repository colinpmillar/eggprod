\name{get_dates}
\alias{get_dates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get dates
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_dates(year = c("1977", "1980", "1983", "1986", "1995", "2001", "2007"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{year}{
which year to use the dates from.
}
}
\details{
  information about survey periods and start and end dates:  
  dates are measured in days since 1 January.
  however, note that start dates correspond to the midnight before, 
  whereas end dates correspond to the midnight after;  
  thus a start date of 10 February is actually midnight between 9 and 10 February and is recorded as 40;  
  an end date of 31 July is actually midnight between 31 July and 1 August and is recorded as 212.
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
function (year = c("1977", "1980", "1983", "1986", "1995", "2001", 
    "2007")) 
{
    year <- match.arg(year)
    if (year == "2007") {
        out <- list(nperiod = 5L, period = cbind(2:6, c(66, 99, 
            127, 155, 176), c(98, 126, 154, 175, 197)), dates = c(40, 
            212))
    }
    else if (year == "2001") {
        out <- list(nperiod = 5L, period = cbind(2:6, c(41, 98, 
            133, 161, 182), c(97, 132, 160, 181, 212)), dates = c(40, 
            212))
    }
    else if (year == "1995") {
        out <- list(nperiod = 5L, period = cbind(3:7, c(41, 104, 
            136, 159, 180), c(103, 135, 158, 179, 212)), dates = c(40, 
            212))
    }
    else if (year == "1986") {
        out <- list(nperiod = 4L, period = cbind(1:4, c(41, 117, 
            155, 182), c(116, 154, 181, 212)), dates = c(40, 
            212))
    }
    else if (year == "1983") {
        out <- list(nperiod = 3L, period = cbind(1:3, c(41, 112, 
            138), c(111, 137, 212)), dates = c(40, 212))
    }
    else if (year == "1980") {
        out <- list(nperiod = 5L, period = cbind(1:5, c(41, 78, 
            105, 135, 187), c(79, 104, 134, 186, 212)), dates = c(40, 
            212))
    }
    else if (year == "1977") {
        out <- list(nperiod = 5L, period = cbind(1:5, c(41, 87, 
            107, 141, 166), c(86, 106, 140, 165, 212)), dates = c(40, 
            212))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }

