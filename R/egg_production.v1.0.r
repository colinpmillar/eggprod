
eggprod <-
function(dates, year, survey_area, what = c("mackerel", "horse mackerel"))
{
  
  what <- match.arg(what)
  nmin <- survey_area $ nmin
  nmax <- survey_area $ nmax
  longmin <- survey_area $ longmin
  longmax <- survey_area $ longmax
  
  names(longmin) <- names(longmax) <- paste(nmin:nmax)  
  
  wmin <- round(2 * min(longmin) + 1)
  wmax <- round(2 * max(longmax))
  
  in_survey_area <- array(FALSE, c(wmax-wmin+1, nmax-nmin+1), list(wmin:wmax, nmin:nmax))
  area <- array(numeric(1), c(wmax-wmin+1, nmax-nmin+1), list(wmin:wmax, nmin:nmax))

  # create a logical array showing whether each rectangle is in the survey area, and calculates the area of that rectangle  
  for (n in nmin:nmax)
  {
    cn <- paste(n)
    if (round(longmax[cn] - longmin[cn]) > 0) # checks if any rectangles present at that latitude
    {
      iw1 <- iwest(longmin[cn] + 0.01)
      iw2 <- iwest(longmax[cn] - 0.01)
      in_survey_area[paste(iw1:iw2), cn] <- TRUE 
      area[paste(iw1:iw2), cn] <- statarea(n)
    }
  }
    
  # calculate daily egg production for each period in turn
  period.data <- 
    lapply(seq.int(dates $ nperiod), 
      function(i)
      {
        # calculate number of hauls and mean number of eggs per rectangle
        # calculates mean number of eggs per rectangle during period
        nhaul <- array(integer(1), c(wmax-wmin+1, nmax-nmin+1), list(wmin:wmax, nmin:nmax))
        mean_eggs <- array(numeric(1), c(wmax-wmin+1, nmax-nmin+1), list(wmin:wmax, nmin:nmax))
                                        
        inputfile <- paste(year, "p", dates$period[i], ".dat", sep = "")   # datafiles of form 1998p1.dat etc
        
        dat <- read.table(inputfile, skip = 1)
        names(dat) <- c("m_eggs", "h_eggs", "lat", "lon")
        dat $ eggs <- if (what == "mackerel") dat $ m_eggs else dat $ h_eggs
        dat $ w <- iwest(dat $ lon)
        dat $ n <- inorth(dat $ lat)
        dat $ period <- i
       
        test <- with(dat, n >= nmin & n <= nmax & w >= wmin & w <= wmax) 
        dat <- subset(dat, test)
        test <- diag(in_survey_area[paste(dat $ w), paste(dat $ n)])
        dat <- subset(dat, test)
        cw <- paste(dat $ w)
        cn <- paste(dat $ n)
              
        for (j in seq_along(dat[[1]]))
        {
          nhaul[cw[j], cn[j]] <- nhaul[cw[j], cn[j]] + 1
          mean_eggs[cw[j], cn[j]] <- mean_eggs[cw[j], cn[j]] + dat $ eggs[j]
        }
       
        area_fill <- afill (in_survey_area, nhaul, area, nmin, nmax, wmin, wmax)
        
        list(nhaul = nhaul, mean_eggs = mean_eggs, area_fill = area_fill, dat = dat)
      })
  
  all.data <- do.call(rbind, lapply(period.data, "[[", "dat"))
  all.data $ id <- with(all.data, paste(floor(w/2), n, period, sep="."))
  tab <- with(subset(all.data, eggs > 0), table(id))
  cv.data <- subset(all.data, eggs > 0 & id %in% names(tab)[tab>1])
  cv.mod <- lm(log(eggs) ~ factor(id), cv.data)
  se <- summary(cv.mod) $ sigma
  cv <- sqrt(exp(se^2) - 1)
      
  eprod.out <-
    sapply(period.data,
      function(x)
      {  
        out <- subset( with(x, 
                         data.frame(nh = c(nhaul), m = c(mean_eggs/nhaul), 
                                    a = c(area), af = c(area_fill))), nh > 0)
        out $ eprod <- with(out, m * a)
        out $ eprodf <- with(out, m * af)
        out $ eprod.var <- with(out, eprod^2/nh)
        out $ eprodf.var <- with(out, eprodf^2/nh)
        
        colSums(out[c("eprod", "eprodf", "eprod.var", "eprodf.var")])
      })
  
  out <- list(what = what, dates = dates)
  out $ eprod <- eprod.out[1,]
  out $ pvar <- eprod.out[3,]
  out $ filleprod <- eprod.out[2,]
  out $ fillpvar <- eprod.out[4,]

  # calculate final production estimate   
  out <- c(with(out, prodcalc (filleprod, fillpvar, dates) ), out)
  
  out $ cv.approx <- se
  out $ cv <- cv
  out $ nzerosq <- length(unique(cv.data $ id))
  out $ nzerodat <- nrow(cv.data)
   
  class(out) <- "eprod"
  out
}

print.eprod <- 
function(x, approx.cv = TRUE, ...)
{
  # write out results
  
  cat(sprintf('Daily egg production estimates of %s\n\n', x $ what))
  cat(sprintf('Period                   Production                  Variance (excluding cv)\n'))
  cat(sprintf('          Raw estimate     Fill-ins        Total       Raw estimate    Total\n'))
  
  for (p in seq.int(x $ dates $ nperiod))
    cat(sprintf('%4i  %14.5e %14.5e %14.5e %14.5e %14.5e\n', 
                x$dates$period[p,1], x$eprod[p], x$filleprod[p] - x$eprod[p], x$filleprod[p], x$pvar[p], x$fillpvar[p]))     
  cat("\n")

  cat(sprintf('period egg production estimates\n\n'))
  pprod <- x $ pprod
  rownames(pprod) <- c("*", c(rbind(paste(x $ dates $ period[,1]), "*")))
  print(pprod)
  cat("\n")
  
  cat(sprintf("Total annual egg production %7.4e\n", x $ taep))
  if (approx.cv)
  {
    cat(sprintf("Variance                    %7.4e\n", with(x, taep.var * cv.approx^2)))
    cat(sprintf("Standard error              %7.4e\n", with(x, sqrt(taep.var) * cv.approx)))
    cat(sprintf("Coefficient of variation    %3.2f\n", with(x, sqrt(taep.var) * cv.approx / taep * 100)))
  } else
  {
    cat(sprintf("Variance                    %7.4e\n", with(x, taep.var * cv^2)))
    cat(sprintf("Standard error              %7.4e\n", with(x, sqrt(taep.var) * cv)))
    cat(sprintf("Coefficient of variation    %3.2f\n", with(x, sqrt(taep.var) * cv / taep * 100)))
  }
  cat(sprintf("\ncv estimation based on %i non zero data points from %i squares", x $ nzerodat, x $ nzerosq))
  cat("\n")
}

# to calculate the half-rectangle west containing longitude
iwest <- function(longitude) floor(2 * longitude) + 1

# to calculate the rectangle north containing latitude
inorth <- function(latitude) floor(2 * latitude) - 71

# calculates the area of a statistical rectangle at latitude n
statarea <- function(n) cos((n + 71.5) * pi / 360.0) * 30^2 * 1853.2^2

afill <-
function (in_survey_area, nhaul, area, nmin, nmax, wmin, wmax)
{
  # calculates the values by which each mean is multiplied, having filled-in

  fill <- ifelse(nhaul == 0, 0, area)
  fill <- cbind(0, rbind(0, fill, 0), 0)
  nhaul2 <- cbind(0, rbind(0, nhaul, 0), 0)
  
  mask <- rbind(c(-1, 1, 0, 0,-1, 1,-1, 1) + 1, 
                c( 0, 0,-1, 1,-1,-1, 1, 1) + 1)

  for (w in 1:(wmax - wmin))
  for (n in 1:(nmax - nmin))
  {
    if (nhaul[w, n] == 0 && in_survey_area[w, n])
    {
      mask2 <- t(mask + c(w, n))
      # calculate number of adjacent, non-diagonal, rectangles with non-zero hauls  
      nct <- sum( pmin(1, nhaul2[ mask2[1:4,] ]) )
      if (nct >= 2)
      {
        # calculate total number of adjacent rectangles with non-zero hauls
        nct <- nct + sum( pmin(1, nhaul2[mask2[5:8,]  ]) )
        # and now fill in those rectangles
        fill[mask2] <- fill[mask2] + ifelse(nhaul2[mask2] > 0, area[w, n] / nct, 0)
      }
    }
  }

  fill[-c(1,nrow(fill)),-c(1,ncol(fill))]
}


prodcalc <-
function (daily_prod, var, dates)
{  
 
  days <- with(dates, diff(c(dates[1], t(period[,-1]), dates[2])))
  nperiod <- dates $ nperiod
  dmat <- c(0, days, 0)[rep(1:5, nperiod) + rep(1:nperiod *2 - 2, each = 5)]
  d <- t(matrix(dmat, 5, nperiod))
  lambda <- d[,3] + d[,2] * (d[,1] + d[,2])/(d[,1] + 2*d[,2] + d[,3]) + 
                    d[,4] * (d[,4] + d[,5])/(d[,3] + 2*d[,4] + d[,5])
  
  sum(lambda * daily_prod)
  
  
  
  np2 <- 2 * nperiod
  
  sampled <- c(FALSE, rep(c(TRUE, FALSE), nperiod)) 

  new_daily_prod <- double(length(days))
  new_daily_prod[sampled] <- daily_prod  # daily production in each sampled period
  var_wt <- days[sampled]

  wk1 <- days[1] / (2 * days[1] + days[2])                                       # daily production between start date and first
  new_daily_prod[1] <- daily_prod[1] * wk1                                       # sampled period
  var_wt[1] <- var_wt[1] + wk1 * days[1]

  for (i in seq(3, np2 - 1, by = 2))
  {                                                           # daily production between sampled periods
    j <- i / 2
    wk1 <- days[i] + days[i+1]
    wk2 <- days[i-1] + days[i]
    new_daily_prod[i] <- (daily_prod[j] * wk1 + daily_prod[j + 1] * wk2) / (wk1 + wk2)
    var_wt[j] <- var_wt[j] + days[i] * wk1 / (wk1 + wk2)
    var_wt[j + 1] <- var_wt[j + 1] + days[i] * wk2 / (wk1 + wk2)
  }                                   

  wk1 <- days[np2+1] / (2 * days[np2+1] + days[np2])                             # daily production between last sampled period and
  new_daily_prod[np2+1] <- daily_prod[nperiod] * wk1                             # end date
  var_wt[nperiod] <- var_wt[nperiod] + wk1 * days[nperiod]

  pd <- new_daily_prod * days                                                    # period production
  pdtot <- sum(pd)                                                               # total production
  vartot <- sum(var * var_wt^2 )                                                 # variance 

  list( pprod = cbind(new_daily_prod, days, pd), taep = pdtot, taep.var = vartot )

}






get_dates <- 
function (year = c("1977", "1980", "1983", "1986", "1995", "2001", "2007"))
{
  
  year <- match.arg(year)
  
  if (year == "2007")
  {
    out <-
     list(nperiod = 5L, 
        period  = cbind(2:6, c( 66,  99, 127, 155, 176 ), c( 98, 126, 154, 175, 197 )), 
        dates   = c(40, 212))
  } else
  if (year == "2001")
  {
    out <-
     list(nperiod = 5L, 
        period  = cbind(2:6, c( 41, 98, 133, 161, 182 ), c( 97, 132, 160, 181, 212 )), 
        dates   = c(40, 212))
  } else
  if (year == "1995")
  {
    out <-
     list(nperiod = 5L, 
        period  = cbind(3:7, c( 41, 104, 136, 159, 180 ), c( 103, 135, 158, 179, 212 )), 
        dates   = c(40, 212))
  } else
  if (year == "1986")
  {
    out <-
     list(nperiod = 4L, 
        period  = cbind(1:4, c( 41, 117, 155, 182 ), c( 116, 154, 181, 212 )), 
        dates   = c(40, 212))
  } else
  if (year == "1983")
  {
    out <-
     list(nperiod = 3L, 
        period  = cbind(1:3, c( 41, 112, 138 ), c( 111, 137, 212 )), 
        dates   = c(40, 212))
  } else
    if (year == "1980")
  {
    out <-
     list(nperiod = 5L, 
        period  = cbind(1:5, c( 41, 78, 105, 135, 187 ), c( 79, 104, 134, 186, 212 )), 
        dates   = c(40, 212))
  } else
  if (year == "1977")
  {
    out <-
     list(nperiod = 5L, 
        period  = cbind(1:5, c( 41, 87, 107, 141, 166 ), c( 86, 106, 140, 165, 212 )), 
        dates   = c(40, 212))
  }

  out
}

get_survey_area <-
function (year = c("1998", "1998_extended", "1995"))
{
  year <- match.arg(year)
  
  if (year == "1995")
  {
    out <-
    list( nmin = 17L, 
          nmax = 44L,
          longmin =                                c( 1.0,  1.0,  1.0,  1.0,      
                  1.5,  2.0,  2.5,  4.0,  4.5,  5.0,  6.0,  6.5,  7.0,  7.5,      
                  8.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,  9.0,  8.0,  8.0,     
                  8.0,  8.0,  8.0,  8.0 ),
          longmax =                                c( 3.0,  5.0,  5.0,  6.0,     
                  6.5,  9.5, 11.0, 12.0, 13.0, 14.0, 14.0, 14.0, 14.0, 14.0,     
                 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 13.5, 11.5, 11.0,     
                 10.5, 10.5, 10.5, 10.5  ))
  } else
  if (year == "1998")
  {
    out <-
    list( nmin = 17L, 
          nmax = 45L,
          longmin =                                c( 1.0,  1.0,  1.0,  1.0,     
                  1.5,  2.0,  2.5,  4.0,  4.5,  5.0,  6.0,  6.5,  7.0,  7.5,      
                  8.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,  9.0,  8.0,  8.0,     
                  8.0,  8.0,  8.0,  8.0,  8.0 ),
          longmax =                                c( 3.0,  5.0,  5.0,  6.5,     
                 10.0, 10.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 16.0,     
                 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 13.5, 11.5, 11.0,     
                 10.5, 10.5, 11.0, 11.0, 11.0  ))
  } else
  if (year == "1998_extended")
  {
    out <-
    list( nmin = 17L, 
          nmax = 46L,
          longmin =                                c( 1.0,  1.0,  1.0,  1.0,       
                  1.5,  2.0,  2.5,  4.0,  4.5,  5.0,  6.0,  6.5,  7.0,  7.5,     
                  8.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,  9.0,  8.0,  8.0,     
                  8.0,  8.0,  8.0,  8.0,  8.0,  8.0 ),
          longmax =                                c( 3.0,  5.0, 20.0, 20.0,     
                 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0,     
                 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0,     
                 20.0, 20.0, 20.0, 20.0, 20.0, 20.0  ))
  }

  out
}
