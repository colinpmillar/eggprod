\name{eggprod}
\alias{eggprod}
\title{
eggprod
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eggprod(dates, year, survey_area)
}
\arguments{
  \item{dates}{
%%     ~~Describe \code{dates} here~~
}
  \item{year}{
%%     ~~Describe \code{year} here~~
}
  \item{survey_area}{
%%     ~~Describe \code{survey_area} here~~
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
function (dates, year, survey_area) 
{
    nmin <- survey_area$nmin
    nmax <- survey_area$nmax
    longmin <- survey_area$longmin
    logmax <- survey_area$longmax
    names(longmin) <- names(longmax) <- paste(nmin:nmax)
    wmin <- round(2 * min(longmin) + 1)
    wmax <- round(2 * max(longmax))
    tmp <- setup(longmin, longmax, nmin, nmax, wmin, wmax)
    in_survey_area <- tmp[[1]]
    area <- tmp[[2]]/1e+10
    eprod.list <- list()
    for (p in seq.int(dates$nperiod)) {
        tmp <- estimate(dates$period[p, 1], year, in_survey_area, 
            nmin, nmax, wmin, wmax)
        nhaul <- tmp[[1]]
        mean_eggs_m <- tmp[[2]]
        mean_eggs_h <- tmp[[3]]
        eprod <- rbind(atotal(nhaul, mean_eggs_m, area), atotal(nhaul, 
            mean_eggs_h, area))
        afillin <- afill(in_survey_area, nhaul, area, nmin, nmax, 
            wmin, wmax)
        eprod <- rbind(eprod, atotal(nhaul, mean_eggs_m, afillin), 
            atotal(nhaul, mean_eggs_h, afillin))
        dim(eprod) <- c(2, 2, 2)
        eprod.list[[p]] <- eprod[1, , ]
    }
    eprod <- sapply(eprod.list, function(x) x[1, 1])
    filleprod <- sapply(eprod.list, function(x) x[2, 1])
    pvar <- sapply(eprod.list, function(x) x[1, 2])
    fillpvar <- sapply(eprod.list, function(x) x[2, 2])
    cat(sprintf("Daily egg production estimates\n\n"))
    cat(sprintf("Period                   Production                  Variance (excluding cv)\n"))
    cat(sprintf("          Raw estimate     Fill-ins        Total       Raw estimate    Total\n"))
    for (p in seq.int(dates$nperiod)) cat(sprintf("\%4i  \%14.5f \%14.5f \%14.5f \%14.5f \%14.5f\n", 
        dates$period[p, 1], eprod[p], filleprod[p] - eprod[p], 
        filleprod[p], pvar[p], fillpvar[p]))
    cat("\n")
    out <- prodcalc(filleprod, fillpvar, dates)
    fname <- ""
    for (i in 1:nrow(out[[1]])) cat(file = fname, sprintf("\%14.3f \%14.1f \%14.3f\n", 
        out[[1]][i, 1], out[[1]][i, 2], out[[1]][i, 3]))
    cat(file = fname, "\n")
    cat(file = fname, sprintf("TAEP NA \%14.3f\n", out[[2]][1]))
    cat(file = fname, sprintf("VAR NA \%14.3f\n", out[[2]][2]))
    cat(file = fname, sprintf("STD NA \%14.3f\n", out[[2]][3]))
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ model }
