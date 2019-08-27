#' distrib
#'
#' produces random samples from distributuions
#'
#' @param shape Required with no default. The permited inputs are Bernoulli, binomial, beta, discrete, empirical,
#' exponential, gamma, lognormal, normal, point, probability, triangle, uniform, and Weibull.
#'
#' @param par1 optional value requared for some shapes. Defulat = none
#'
#' @param par2 optional value requared for some shapes. Defulat = none
#'
#' @param par3 optional value requared for some shapes. Defulat = none
#'
#' @param par4 optional value requared for some shapes. Defulat = none
#'
#' @param It optional value requared for some shapes. Defulat = none
#'
#' @param Ut optional value requared for some shapes. Defulat = none
#'
#' @param resamp optional value requared for some shapes. Defulat = 'y'
#'
#' @param n Optional	Default = 1
#'
#' @param q	Optional	Default = none
#'
#' @param p	Optional	Default = c(1)
#'
#' @param  v Optional	Default = none
#'
#' @details  This process is central to SHEDS because it is a stochastic model.  The "shape" is the essential argument,
#' with others being required for certain shapes.  For all shapes, one of  "n" or "q" must be specified.  If "n" is given,
#' then Distrib returns a vector of "n" independent random samples from the specified distribution.  If "q" is given, it
#' must be a vector of numeric values, each between zero and one.  These are interpreted as the quantiles of the distribution
#' to be returned.  When "q" is given, Distrib does not generate any random values, it just evaluates the requested quantiles.
#' The empirical shape requires argument "v" as a list of possible values to be returned, each with equal probability.
#' The other shapes require one or more of par1-par4 to be specified.  See the SHEDS Technical Manual for more details on the
#' meanings of par1-par4, which vary by shape.  "Lower.trun" is the lower truncation point, meaning the smallest value that can
#' be returned.  Similarly, "upper.trun" is the largest value that may be returned.  Not all distributions use lower.trun and/or
#' upper.trun, but they should be specified for unbounded shapes like the Normal distribution.  "Resamp" is a flag to indicate
#' the resolution for generated values outside the truncation limits.  If resamp="yes" then effectively new values are generated
#' until they are within the limits.  If resamp="no", values outside the limits are moved to those limits.  "P" is a list of
#' probabilities that are used only with the "discrete" or "probability" shapes.  The "p" values are essentially weights for a
#' list of discrete values that may be returned.  The "empirical" distribution also returns discrete values, but assigns them
#' equal weights, so then "p" is not needed.
#'
#' @return A vector of "n" values from one distribution, where "n" is
#' either the input argument (if given), or the length of the input vector "q".
#'
#' @author  Kristin Isaacs, Graham Glen
#'
#' @export

distrib = function(shape="",par1=NA,par2=NA,par3=NA,par4=NA,lt=NA,
                   ut=NA, resamp="y",n=1,q=NA,p=c(1),v="" ) {
  # Distrib generates samples from a distribution. If a vector of quantiles
  # q is given, the corresponding values are returned. Otherwise n determines
  # the number of samples, but then quantiles are randomly generated first.
  # Written by Graham Glen, Alion Science and Technology, July 2012.
  # Last modified by GG (ICF International) on August 28, 2016.

  m <- ""
  n <- round(n)
  if (n<=0) m <- paste("# samples requested = ",n)
  if (length(par1)>1|length(par2)>1|length(par3)>1|length(par4)>1) {
    m <- "Non-scalar distribution parameters"
  }
  if (is.na(q[1])) q <- runif(n)
  s <- strtrim(tolower(shape),4)
  r <- strtrim(tolower(resamp),1)
  if (is.na(lt)&is.na(ut)) r <- "n"
  if (!is.na(ut)&!is.na(lt)&ut<lt) m <- paste("Truncation limits",lt,"and",ut)
  if (!is.numeric(q)) m <- "Quantiles are not numeric"
  if (min(q)<0 || max(q)>1) m <- "Quantiles not all between 0 and 1"

  if (m!="") {
  } else if (s=="bern" || s=="bino") {
    if (is.na(par1)) par1 <- 0.5
    lt <- NA
    ut <- NA
    if (par1<0 | par1>1) { m <- paste("Invalid binomial parameter",par1)
    } else {
      x <- q
      x[q<=1-par1] <- 0
      x[q> 1-par1] <- 1
    }
  } else if (s=="beta") {
    if (is.na(par1)) par1 <- 1
    if (is.na(par2)) par2 <- 1
    if (is.na(par3)&&is.na(par4)) { par3 <- 0; par4 <- 1
    } else if (is.na(par4)) { par4 <- par3+1
    } else if (is.na(par3)) par3 <- par4-1
    if (!is.na(lt) && lt>par4) m <- "Lower beta truncation above par4"
    if (!is.na(ut) && ut<par3) m <- "Upper beta truncation below par3"
    if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Lower truncation above upper"
    w <- par4-par3
    if (par1<=0 | par2<=0 | w<0) {
      m <- paste("Invalid beta parameters",par1,par2,par3,par4)
    } else if (r=="y") {
      if (!is.na(lt)) qlo <- pbeta((lt-par3)/w,par1,par2) else qlo <- 0
      if (!is.na(ut)) qhi <- pbeta((ut-par3)/w,par1,par2) else qhi <- 1
      q <- qlo+q*(qhi-qlo)
    }
    if (m=="") x <- par3 + w*qbeta(q,par1,par2)
  } else if (s=="disc" || s=="prob") {
    if (mode(p)=="character") p <- scan(text=p,quiet=TRUE)
    if (mode(v)!="numeric") {
      if (mode(try(scan(text=v,quiet=TRUE),silent=TRUE))=="numeric") {
        v <- c(scan(text=v,quiet=TRUE),recursive=TRUE)
      } else {
        v <- c(scan(what=list(""),text=v,quiet=TRUE),recursive=TRUE)
      }
    }
    if (length(v)==0) v <- 0
    if (s=="prob") v <- 1:length(p)
    if (s=="disc" && length(p)==1 && length(v)>1) p <- rep(1,length(v))
    if (s=="disc" && length(v)==1 && length(p)>1) v <- 1:length(p)
    if (!is.na(lt)) lt <- min(v[which(v>=lt)])
    if (!is.na(ut)) ut <- max(v[which(v<=ut)])
    if (length(p)!=length(v)) {
      m <- paste("Unequal vectors, v=",length(v),", p=",length(p))
    } else if (min(p)<0) { m <- "Negative probabilities found"
    } else if (r=="y") {
      if (!is.na(lt)) p[v<lt] <- 0
      if (!is.na(ut)) p[v>ut] <- 0
    }
    if (length(p[p>0])==0) m <- "All probabilities are zero"
    if (m=="") {
      t <- cumsum(p/sum(p))
      t[t>0.99999999] <- 1.00000001
      x <- v[mapply(function(q,t) {which.max(cummax(q<t))},q,MoreArgs=list(t))]
    }
  } else if (s=="empi") {
    if (mode(v)!="numeric") {
      if (mode(try(scan(text=v,quiet=TRUE),silent=TRUE))=="numeric") {
        v <- c(scan(text=v,quiet=TRUE),recursive=TRUE)
      } else {
        v <- c(scan(what=list(""),sep=",",text=v,quiet=TRUE),recursive=TRUE)
      }
      lt<-NA
      ut<-NA
    }
    if (length(v)==0) v <- 0
    if (!is.na(lt)) lt <- min(v[which(v>=lt)])
    if (!is.na(ut)) ut <- max(v[which(v<=ut)])
    if (r=="y") {
      if (!is.na(lt)) v <- v[v>=lt]
      if (!is.na(ut)) v <- v[v<=ut]
    }
    if (length(v)==0) m <- "No empirical values"
    if (m=="") x <- v[round(0.5+length(v)*q)]
  } else if (s=="expo") {
    if (is.na(par1)) par1 <- 1
    if (is.na(par2)) par2 <- 0
    if (!is.na(lt) && lt<par2) lt <- par2
    if (!is.na(ut) && ut<par2) m <- "Upper expo. truncation below par2"
    if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Lower truncation above upper"
    if (par1<=0) { m <- paste("Invalid exponential parameter",par1)
    } else if (r=='y') {
      if (!is.na(lt)) qlo <- pexp(lt-par2,par1)  else qlo <- 0
      if (!is.na(ut)) qhi <- pexp(ut-par2,par1)  else qhi <- 1
      q <- qlo+q*(qhi-qlo)
    }
    if (m=="") x <- par2 + qexp(q,par1)
  } else if (s=="gamm") {
    if (is.na(par1)) par1 <- 1
    if (is.na(par2)) par2 <- 1
    if (is.na(par3)) par3 <- 0
    if (!is.na(lt) && lt<par3) lt <- par3
    if (!is.na(ut) && ut<par3) m <- "Upper gamma truncation below par3"
    if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Upper truncation above lower"
    if (par1<=0 | par2<=0) {
      m <- paste("Invalid gamma parameters",par1,par2,par3)
    } else if (r=='y') {
      if (!is.na(lt)) qlo <- pgamma(lt-par3,par1,1/par2,par2) else qlo <- 0
      if (!is.na(ut)) qhi <- pgamma(ut-par3,par1,1/par2,par2) else qhi <- 1
      q <- qlo+q*(qhi-qlo)
    }
    if (m=="") x <- par3 + qgamma(q,par1,1/par2,par2)
  } else if (s=="logn") {
    if (is.na(par1)) par1 <- 1
    if (is.na(par2)) par2 <- exp(1)
    if (is.na(par3)) par3 <- 0
    if (!is.na(lt) && lt<par3) lt <- par3
    if (!is.na(ut) && ut<par3) m <- "Upper lognormal truncation below par3"
    if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Upper truncation above lower"
    if (par1<=0 | par2<1) {
      m <- paste("Invalid lognormal parameters",par1,par2,par3)
    } else if (r=='y' && par2>1) {
      if (!is.na(lt)) qlo <- plnorm(lt-par3,log(par1),log(par2)) else qlo <- 0
      if (!is.na(ut)) qhi <- plnorm(ut-par3,log(par1),log(par2)) else qhi <- 1
      q <- qlo+q*(qhi-qlo)
    }
    if (m=="") x <- par3 + qlnorm(q,log(par1),log(par2))
  } else if (s=="norm") {
    if (is.na(par1)) par1 <- 0
    if (is.na(par2)) par2 <- 1
    if (par2<0) {m <- paste("Invalid normal parameters",par1,par2)
    if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Lower truncation above upper"
    } else if (r=='y' && par2>0) {
      if (!is.na(lt)) qlo <- pnorm(lt,par1,par2) else qlo <- 0
      if (!is.na(ut)) qhi <- pnorm(ut,par1,par2) else qhi <- 1
      q <- qlo+q*(qhi-qlo)
    }
    if (m=="") x <- qnorm(q,par1,par2)
  } else if (s=="poin") {
    if (is.na(par1)) par1 <- 0
    lt <- NA
    ut <- NA
    x <- rep(par1,length(q))
  } else if (s=="tria") {
    if (is.na(par1)&&is.na(par2)) { par1 <- 0; par2 <- 1
    } else if (is.na(par2)) { par2 <- par1+1
    } else if (is.na(par1))   par1 <- par2-1
    if (is.na(par3)) par3 <- (par1+par2)/2
    if (par3>par2) { t<-par2; par2<-par3; par3<-t }
    if (!is.na(lt) && lt>par2) m <- "Lower triangle truncation above par2"
    if (!is.na(ut) && ut<par1) m <- "Upper triangle truncation below par1"
    if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Upper truncation above lower"
    if (par1>par2 | par3<par1 | par3>par2 ) {
      m <- paste("Invalid triangle parameters",par1,par2,par3)
    }
    if (par1==par2) return(rep(par1,length(q)))
    p <- (par3-par1)/(par2-par1)
    if (r=='y') {
      if (!is.na(lt) && lt>par1 && lt<=par2) {
        if (lt==par3) { qlo <- p
          } else if (lt<par3) { qlo <-   (lt-par1)^2/((par2-par1)*(par3-par1))
          } else if (lt>par3) { qlo <- 1-(par2-lt)^2/((par2-par1)*(par2-par3))
        }
      } else qlo <- 0
      if (!is.na(ut) && ut>=par1 && ut<par2) {
        if (ut==par3) {qhi <- p
          } else if (ut<par3) { qhi <-   (ut-par1)^2/((par2-par1)*(par3-par1))
          } else if (ut>par3) { qhi <- 1-(par2-ut)^2/((par2-par1)*(par2-par3))
        }
      } else qhi <- 1
      q <- qlo+q*(qhi-qlo)
    }
    if (m=="") {
      x  <- par1 + sqrt(   q *(par2-par1)*(par3-par1))
      x2 <- par2 - sqrt((1-q)*(par2-par1)*(par2-par3))
      x[x2>par3] <- x2[x2>par3]
    }
  } else if (s=="unif") {
    if (is.na(par1)) par1 <- 0
    if (is.na(par2)) par2 <- 1
    if (par2<par1) m <- paste("Invalid uniform parameters", par1,par2)
    if (!is.na(lt) && lt>par2) m <- "Lower uniform truncation above par2"
    if (!is.na(ut) && ut<par1) m <- "Upper uniform truncation below par1"
    if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Upper truncation above lower"
    if (r=='y' && par2>par1) {
      if (!is.na(lt)) qlo <- punif(lt,par1,par2) else qlo <- 0
      if (!is.na(ut)) qhi <- punif(ut,par1,par2) else qhi <- 1
      q <- qlo+q*(qhi-qlo)
    }
    if (m=="") x <- qunif(q,par1,par2)
  } else if (s=="weib") {
    if (is.na(par1)) par1 <- 1
    if (is.na(par2)) par2 <- 1
    if (is.na(par3)) par3 <- 0
    if (!is.na(lt) && lt<par3) lt <- par3
    if (!is.na(ut) && ut<par3) m <- "Upper Weibull truncation below par3"
    if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Upper truncation above lower"
    if (par1<=0 | par2<=0) {
      m <- paste("Invalid Weibull parameters",par1,par2,par3)
    }
    if (r=='y') {
      if (!is.na(lt)) qlo <- pweibull(lt-par3,par1,par2) else qlo <- 0
      if (!is.na(ut)) qhi <- pweibull(ut-par3,par1,par2) else qhi <- 1
      q <- qlo+q*(qhi-qlo)
    }
    if (m=="") x <- par3+qweibull(q,par1,par2)
  } else m <- paste("Unknown Distrib shape ",s)

  if (m != "") {cat("\nError in distrib: ",m,"\n"); return(NULL)
  } else {
    if (!is.na(lt)) x <- mapply(max,lt,x)
    if (!is.na(ut)) x <- mapply(min,ut,x)
    return(x)
  }
}

#' p.round
#'
#' p.round performs stochastic or probabilistic rounding of non-integer values.
#'
#' @param x	Required 	Default = none
#'
#' @param q	Optional	Default = none
#'
#' @details The input "x" is a vector of values to be rounded.  Each value of x is rounded independently, either up or down.
#' The fractional value of x is used to determine weights for rounding in each direction.  For example, if x=4.3, then it is
#' rounded up to 5 with probability 0.3, else rounded down to 4 (with probability 0.7).   The p.round function always returns
#' integer values, and does not introduce bias.  In the above example, if the argument were 4.3 a large number of times, then
#' the returned values (all either 4 or 5) would average 4.3 with a small residual error which approaches zero as "n" gets large.
#'
#' @return A vector of the same length as the input "x", containing all integer values.
#'
#' @author  Kristin Isaacs, Graham Glen
#'
#' @export

p.round = function(x=NULL,q=NA) {
  # A probabilistic rounding function, by Graham Glen, Sep 2013
  if (is.na(q[1])) q <- runif(length(x))
  y   <- floor(x)
  z   <- y + as.numeric(q > 1+y-x)
  return(z)
}

#' quantiles
#'
#' This function is similar to the built-in R function "quantile", but it returns a list of pre-selected quantiles of the set
#' of values in the vector "x".
#'
#' @param x	Default = none
#'
#' @details This function is similar to the built-in R function "quantile", but it returns a list of pre-selected quantiles
#' of the set of values in the vector "x".  Specifically, it returns all the following quantiles:
#' .005, .01, .025, .05, .1, .15, .2, .25, .3, .4, .5, .6, .7, .75, .8, .85, .9, .95, .975, .99, .995
#'
#' @return This function is used to construct the tables in the "Allstats" output files.
#'
#' @author  Kristin Isaacs, Graham Glen
#'
#' @export

quantiles = function(x) {
  v <- c(.005,.01,.025,.05,.1,.15,.2,.25,.3,.4,.5,.6,
         .7,.75,.8,.85,.9,.95,.975,.99,.995)
  y <- quantile(x,v,type=3,na.rm=TRUE)
  return(y)
}

#' summarize.chemcial
#'
#' Summarize.chemical writes a .csv file containing a summary of the exposure and dose results from the object "x".
#'
#' @param x	Default = none		Exposure data set
#'
#' @param c	Default = none	Index # for chemical
#'
#' @param chem Default = none	CAS for chemical
#'
#' @param chemical Default = none	Full chemical name
#'
#' @param set	Default = none	Index # for set of simulated persons
#'
#' @param sets Default = none		Total # of sets in this SHEDS run
#'
#' @param specs		Default = none List of settings from the "run" input file
#'
#' @details Tables are produced for each of the following cohorts: males, females, females ages 16-49, age 0-5, age 6-11,
#' age 12-19, age 20-65, age 66+, and a table for all persons. The variables that are summarized are: dermal exposure,
#' ingestion exposure, inhalation exposure, inhaled dose, intake dose, dermal absorption, ingestion absorption, inhalation
#' absorption, total absorption in micrograms per day, total absorption in milligrams per kilogram per day, and chemical mass
#' down the drain.
#' In each table, the following statistics are computed across the appropriate subpopulation: mean, standard deviation, and
#' quantiles .005, .01, .025, .05, .1, .15, .2, .25, .3, .4, .5, .6, .7, .75, .8, .85, .9, .95, .975, .99, .995.
#'
#' @return  If "x" is a single set of data (that is, if the argument "set" is between 1 and sets, inclusive), then the .csv
#' file created by this function has the suffix "_set#stats.csv", where "#" is the set number.  If the "set" argument is
#' "allstats", then the .csv file has the suffix _allstats.csv".  The output file contains the tables for all cohorts with a
#' non-zero population.
#'
#' @author  Kristin Isaacs, Graham Glen
#'
#' @export

summarize.chemical = function(x,c,chem,chemical,set,sets,specs) {
  suffix <- paste0("_set",set,"stats")
  mode(x$age) <- "numeric"
 if (set=="allstats") { names(x) <- c("person","gender","age","season","weekend","weight",
     "exp.dermal.tot","exp.inhal.tot","exp.ingest.tot","exp.diet","exp.nondiet","exp.ddd.tot",
     "dose.inhal.tot","dose.intake.ug","dose.intake.mgkg","abs.dermal.ug","abs.inhal.ug",
     "abs.ingest.ug","abs.tot.ug","abs.tot.mgkg","urine.tot.ug","exp.window","conc.inhal.max.prod.aer","conc.inhal.max.prod.vap","exp.inhal.indir","abs.hm.ug")
     suffix <- "_allstats"
  }
  Statistic <- as.list(cbind("Cohort","ug/day","ug/day","ug/m3","ug/day","mg/kg/day","ug/day",
                             "ug/day","ug/day","ug/day","mg/kg/day","g/day","ug/m3","ug/m3","ug/m3","ug/day"))
  a <- rbind(Statistic,cbind("Total",summary.stats(x)))
  b  <- as.data.table(a)
  y <- x[x$gender=="M"]
  if (nrow(y)>0) a <- rbind(a,cbind("Males",summary.stats(y)))
  y <- x[x$gender=="F"]
  if (nrow(y)>0) a <- rbind(a,cbind("Females",summary.stats(y)))
  y <- x[x$gender=="F" & x$age>= 16 & x$age<=49]
  if (nrow(y)>0) a <- rbind(a,cbind("Females_Repro",summary.stats(y)))
  y <- x[x$age<6]
  if (nrow(y)>0) a <- rbind(a,cbind("Age_0_5",summary.stats(y)))
  y <- x[x$age>=6 & x$age<=11]
  if (nrow(y)>0) a <- rbind(a,cbind("Age_6_11",summary.stats(y)))
  y <- x[x$age>=12 & x$age<=19]
  if (nrow(y)>0) a <- rbind(a,cbind("Age_12_19",summary.stats(y)))
  y <- x[x$age>=20 & x$age<=65]
  if (nrow(y)>0) a <- rbind(a,cbind("Age_20_65",summary.stats(y)))
  y <- x[x$age>=66]
  if (nrow(y)>0) a <- rbind(a,cbind("Age_66+",summary.stats(y)))

  if(specs$details!="0" & tolower(specs$details)!="no") {
    n1 <- nrow(b)-1
    stat  <- "pop.mean"
    exder <- signif(as.numeric(b$exp.dermal[n1]),6)
    exing <- signif(as.numeric(b$exp.ingest[n1]),6)
    exinh <- signif(as.numeric(b$exp.inhal[n1]),6)
    doinh <- signif(as.numeric(b$dose.inhal[n1]),6)
    intak <- signif(as.numeric(b$dose.intake[n1]),6)
    abder <- signif(as.numeric(b$abs.dermal.ug[n1]),6)
    abing <- signif(as.numeric(b$abs.ingest.ug[n1]),6)
    abinh <- signif(as.numeric(b$abs.inhal.ug[n1]),6)
    atot1 <- signif(as.numeric(b$abs.tot.ug[n1]),6)
    atot2 <- signif(as.numeric(b$abs.tot.mgkg[n1]),6)
    ddd   <- signif(as.numeric(b$ddd.mass[n1]),6)
    out <- data.frame(stat,exder,exing,exinh,doinh,intak,abder,abing,abinh,atot1,atot2,ddd)
    setnames(out,names(out),c("stat","exp.der","exp.ing","exp.inh","dos.inh","dos.intak","abs.der","abs.ing","abs.inh","tot.ug","tot.mgkg","ddd.mass"))
    cat("\n set=",set,"/",sets," chem=",c,"/",specs$n.chem,"   ",chem,"  ",chemical,"\n")
    print(out)
  }

  dir  <- paste0("output/",specs$run.name)
  name <- paste0(dir,"/CAS_",chem,suffix,".csv")

  write.csv(rbind(a,Statistic),name)
}

#' summary.stats
#'
#' Summary.stats constructs the table of exposure and dose statistics for cohort entered into \code{\link{summarize.chemical}}.
#'
#' @param x. Data set passed from summarize.chemical for a cohort
#'
#' @details Summary.stats is called by \code{\link{summarize.chemical}}.  The input data set "y"
#' is one population cohort from the exposure data set passed into \code{\link{summarize.chemical}}.
#'
#' @return y A data frame object with 23 rows and 11 columns, with each column being an exposure or dose variable,
#' and each row containing a statistic for that variable. For each expsoure varianle the total expsoure, quantiles, mean, and
#' SD.
#'
#' @author  Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{summarize.chemical}}
#'
#' @export

summary.stats = function(x) {
  y1                <- as.numeric(x$exp.dermal.tot)
  exp.dermal        <- quantiles(y1)
  n1                <- length(exp.dermal)+1
  n2                <- n1+1
  exp.dermal[n1]    <- mean(y1,na.rm=TRUE)
  exp.dermal[n2]    <- sd(y1,na.rm=TRUE)
  names(exp.dermal)[n1] <- "mean"
  names(exp.dermal)[n2] <- "sd"
  y2                <- as.numeric(x$exp.ingest.tot)
  exp.ingest        <- quantiles(y2)
  exp.ingest[n1]    <- mean(y2,na.rm=TRUE)
  exp.ingest[n2]    <- sd(y2,na.rm=TRUE)
  y3                <- as.numeric(x$exp.inhal.tot)
  exp.inhal         <- quantiles(y3)
  exp.inhal[n1]     <- mean(y3,na.rm=TRUE)
  exp.inhal[n2]     <- sd(y3,na.rm=TRUE)
  y4                <- as.numeric(x$dose.inhal.tot)
  dose.inhal        <- quantiles(y4)
  dose.inhal[n1]    <- mean(y4,na.rm=TRUE)
  dose.inhal[n2]    <- sd(y4,na.rm=TRUE)
  y5                <- as.numeric(x$dose.intake.mgkg)
  dose.intake       <- quantiles(y5)
  dose.intake[n1]   <- mean(y5,na.rm=TRUE)
  dose.intake[n2]   <- sd(y5,na.rm=TRUE)
  y6                <- as.numeric(x$abs.dermal.ug)
  abs.dermal.ug     <- quantiles(y6)
  abs.dermal.ug[n1] <- mean(y6,na.rm=TRUE)
  abs.dermal.ug[n2] <- sd(y6,na.rm=TRUE)
  y7                <- as.numeric(x$abs.ingest.ug)
  abs.ingest.ug     <- quantiles(y7)
  abs.ingest.ug[n1] <- mean(y7,na.rm=TRUE)
  abs.ingest.ug[n2] <- sd(y7,na.rm=TRUE)
  y8                <- as.numeric(x$abs.inhal.ug)
  abs.inhal.ug      <- quantiles(y8)
  abs.inhal.ug[n1]  <- mean(y8,na.rm=TRUE)
  abs.inhal.ug[n2]  <- sd(y8,na.rm=TRUE)
  y9                <- as.numeric(x$abs.tot.ug)
  abs.tot.ug        <- quantiles(y9)
  abs.tot.ug[n1]    <- mean(y9,na.rm=TRUE)
  abs.tot.ug[n2]    <- sd(y9,na.rm=TRUE)
  y10               <- as.numeric(x$abs.tot.mgkg)
  abs.tot.mgkg      <- quantiles(y10)
  abs.tot.mgkg[n1]  <- mean(y10,na.rm=TRUE)
  abs.tot.mgkg[n2]  <- sd(y10,na.rm=TRUE)
  y11               <- as.numeric(x$exp.ddd.tot)
  ddd.mass          <- quantiles(y11)
  ddd.mass[n1]      <- mean(y11,na.rm=TRUE)
  ddd.mass[n2]      <- sd(y11,na.rm=TRUE)
  y12               <- as.numeric(x$conc.inhal.max.prod.aer)
  conc.max.prod.aer          <- quantiles(y12)
  conc.max.prod.aer[n1]      <- mean(y12,na.rm=TRUE)
  conc.max.prod.aer[n2]      <- sd(y12,na.rm=TRUE)
  y13               <- as.numeric(x$conc.inhal.max.prod.vap)
  conc.max.prod.vap          <- quantiles(y13)
  conc.max.prod.vap[n1]      <- mean(y13,na.rm=TRUE)
  conc.max.prod.vap[n2]      <- sd(y13,na.rm=TRUE)
  y14               <- as.numeric(x$exp.inhal.indir)
  exp.inhal.indir          <- quantiles(y14)
  exp.inhal.indir[n1]      <- mean(y14,na.rm=TRUE)
  exp.inhal.indir[n2]      <- sd(y14,na.rm=TRUE)
  y15               <- as.numeric(x$abs.hm.ug)
  abs.hm.ug         <- quantiles(y15)
  abs.hm.ug[n1]      <- mean(y15,na.rm=TRUE)
  abs.hm.ug[n2]      <- sd(y15,na.rm=TRUE)
  return(cbind(exp.dermal,exp.ingest,exp.inhal,
             dose.inhal,dose.intake,abs.dermal.ug,abs.ingest.ug,
             abs.inhal.ug,abs.tot.ug,abs.tot.mgkg,ddd.mass,conc.max.prod.aer,conc.max.prod.vap,exp.inhal.indir,abs.hm.ug))
}

#' Trimzero
#'
#' This function removes initial zeroes from CAS numbers.
#'
#' @param x aCAS number. Defult is none
#'
#' @param y A dummy argument This helps with the removal of the zeros inteh CAS number
#'
#' @details Each CAS number has three parts, separated by underscores.  The first part is up to seven digits, but optionally,
#' leading zeroes are omitted.  For example, formaldehyde may be either "0000050_00_0" or "50_00_0". SHEDS needs to match
#' CAS numbers across input files, and trimzero is used to ensure matching even when the input files follow different conventions.
#'
#' @return y a shorter CAS number. If the initial part is all zero (as in "0000000_12_3"), one zero is left in the first part
#' (that is, "0_12_3" for this example).
#'
#' @author  Kristin Isaacs, Graham Glen
#'
#' @export

trimzero = function(x,y) {
  y <- x
  for (j in 1:length(x)) {
    i <- 1
    while(substr(x[j],i,i)=="0" & substr(x[j],i+1,i+1)!="_") {
      substr(y[j],i,i)<- " "
      i <- i+1
    }
  }
  return(str_trim(y))
}

#' unpack
#'
#' @details This function is used when ShedsHT is run as an R package. Each time a new working directory is chosen, use unpack()
#' to convert the R data objects into CSV files. Note that both /inputs and /output folders are needed under the chosen
#' directory.Unpack() may be used with a list of object names, in which case just those objects are converted to CSV. This is
#' useful when re-loading one or more defaults into a folder where some of the CSV files have changes, and should not be
#' overwritten. A blank argument or empty list means that all the csv and TXT files in the R package are converted.
#'
#' @author  Kristin Isaacs, Graham Glen
#'
#' @export

unpack = function(filelist="") {
  if (filelist=="") {
    filelist <- c("activity_diaries","chem_props","diet_diaries","exp_factors","fugacity","media","physiology","population",
                  "run_artsandcrafts","run_food_residue","run_other_sources","run_products",
                  "source_chem_ac","source_chem_food","source_chem_products","source_chem_others",
                  "source_scen_food","source_scen_products","source_scen_others",
                  "source_vars_products","source_vars_others","README","run_CPDat")
  }
  inlib <- paste0(getwd(),"/inputs")
  if (dir.exists(inlib)) {

  for (i in 1:length(filelist)) {
    infile  <- filelist[[i]]
    extension <- ".csv"
    if (tolower(substr(infile,1,3))=="run") extension <- ".txt"
    if (tolower(substr(infile,1,3))=="rea") extension <- ".txt"
    csvfile <- paste0(tolower(infile),extension)
       name<-system.file("extdata", csvfile, package = "ShedsHT")
       cat("\n Unpacking ",name)
       file.copy(name, inlib,overwrite=T)
       #ExportDataTables(infile,inlib,csvfile)
  }
  }
  if (!dir.exists(inlib)) {
  cat("\n The Inputs directory does not exist in the SHEDS home directory. Unpack failed.")
  }
}
#' vpos
#'
#' This function locates an item in a list.
#'
#' @details this one was written because it does not require additional R packages and its behavior can be easily examined.
#'
#' @author  Kristin Isaacs, Graham Glen
#'
#' @export

vpos = function(v,list) {
  for (i in 1:length(list)) {
     if (list[i]==v) return(i)
  }
  cat ("\n Variable ",v," not found in list \n")
  return(0)
}

#' write.persons
#'
#' This function writes demographic, exposure, and dose variables to a separate output file for each chemical.
#'
#' @details  THis function rounds the variables to a reasonable precision so that the files are more readable, as unrounded
#' output is cluttered with entries with (say) 14 digits, most of which are not significant.
#'
#' @author  Kristin Isaacs, Graham Glen
#'
#' @export

write.persons = function(x,chem,set,specs) {
  dir    <- paste0("output/",specs$run.name)
  name   <- paste0(dir,"/CAS_",chem,"_all.csv")
  person <- x$person + (set-1)*specs$set.size
  a      <- as.data.table(cbind(person,as.character(x$gender),x$age,
                           as.character(x$season), x$weekend, round(x$weight,3),
                           round(x$exp.dermal.tot,6), round(x$exp.inhal.tot,6),
                           round(x$exp.ingest.tot,6), round(x$exp.dietary.tot,6),
                           round(x$exp.nondiet.tot,6), round(x$exp.ddd.tot,6),
                           round(x$dose.inhal.tot,6), round(x$dose.intake.ug,6),
                           round(x$dose.intake.mgkg,9), round(x$abs.dermal.ug,6),
                           round(x$abs.inhal.ug,6), round(x$abs.ingest.ug,6),
                           round(x$abs.tot.ug,6), round(x$abs.tot.mgkg,9),
                           round(x$urine.tot.ug,6), round(x$exp.window.tot,6),
                           round(x$conc.inhal.max.prod.aer,6),round(x$conc.inhal.max.prod.vap,6),round(x$exp.inhal.indir,6),
                           round(x$abs.hm.ug,6)))
  setnames(a,names(a),c("person","gender","age","season","weekend","weight",
                        "exp.dermal","exp.inhal","exp.ingest","exp.diet","exp.nondiet",
                        "exp.drain","dose.inhal.ug","dose.intake.ug","dose.intake.mgkg",
                        "abs.dermal.ug","abs.inhal.ug","abs.ingest.ug","abs.tot.ug",
                        "abs.tot.mgkg","urine.tot.ug","exp.window","conc.inhal.max.prod.aer","conc.inhal.max.prod.vap","exp.inhal.indir","abs.hm.ug"))
  if (set==1) write.csv(a,name,row.names=FALSE) else {
    write.table(a,name,append=TRUE,sep=',',row.names=FALSE,col.names=FALSE)
  }
}

# Notes for SHEDS-HT Utility module
#
# distrib()   The distrib function generates random samples from distributions. It also
#             can transform quantiles into the appropriate values, given the distribution.
#             Each call can return multiple values, but only from one distribution per
#             call. In SHEDS, with N persons run simultaneously in a set, distrib generates
#             N samples from each distribution. In some cases, not all of them are needed.
#
#             If a vector of quantiles is supplied as a calling argument, then distrib
#             locates those quantiles and returns their values.  If the quantile vector
#             is missing (or specfically, if its first element is NA), then a new vector
#             of quantiles is randomly generated. It then proceeds as if the vector had
#             been supplied.
#
#             The string "m" is initially null, and exists to report any messages that
#             are produced when problems are detected. The code allows just one such
#             message, as each overwrites any existing one. This is mainly for debugging,
#             as SHEDS should not produce any problems here if everything is in order.
#
#             The string "s" contains the first four letters of the distribution type.
#             Most of the code consists of a series of IF...ELSE statements selecting
#             for the value of "s", as the process of transforming the quantiles is
#             different for every distribution type.
#
#             The bernoulli and binomial distributions in SHEDS are the same. Each returns
#             a value of 1 with probability par1, or returns zero otherwise,
#
#             For the beta distribution, par1 and par2 are the shape parameters, par3 is the
#             lower bound, and par4 is the upper bound.  If par3 and par4 are both undefined,
#             they are assumed to be 0 and 1. If one of them is defined, the beta is assumed
#             to have a standard width of one.
#
#             There are three types of discrete distributions in SHEDS, called discrete,
#             probability, and empirical.  The discrete type is the most generall; it requires
#             a vector "p" of probabilities and an equally long vector "v" of possible results.
#             If the probabilities do not sum to one, they are standardized by dividing each
#             by the sum. They are converted to cumulative probabilities, and compared to the
#             vector of selected quantiles, and the value corresponding to the location of the
#             first cumulative probability to exceed each selected quantile is returned.
#
#             For the probability type, the returned values are assumed to be counting numbers
#             1 to N, for a vector of N probabilities.  For the empirical type, the values
#             are needed, but the probability vector "p" is not, as all values are assumed
#             to have equal selection probability.
#
#             For the exponential, par1 is the decay rate (the inverse of the scale parameter),
#             and par2 is the shift (assumed to be zero if absent).
#
#             For the gamma, par1 is the shape and par2 is the scale parameter.  Par3 is optional
#             and is the shift, assumed to be zero if not defined. Par4 is not used.
#
#             For the lognormal, par1 is the geometric mean (which is also the median) and par2
#             is the geometric standard deviation.  These are before any shift is made (par3 is
#             the shift and is assumed to be zero, if missing).  Par1 must be positive and par2
#             must be greater than one.
#
#             For the normal, par1 is the mean and par2 is the standard deviation.
#
#             For the point, par1 is the value that is always returned.
#
#             For the triangle distribution, SHEDS accepts either par1=min, par2=max, par3=peak,
#             or else par1=min, par2=peak, par3=max.  The larger of par2 and par3 is automatically
#             taken to be the max.
#
#             For the uniform, par1 is the minimum and par2 is the maximum,
#
#             For the Weibull distribution, par1 is the shape and par2 is the scale parameter.
#             Par3 is the shift parameter which is optional and has a default value of zero.
#
#             For all distributions, if the resample option r=="y" then the quantiles are mapped
#             onto the range qlo to qhi by q = qlo + (qhi-qlo)*q.  If resample=no then the original
#             quantiles are not modified.  At the end of distrib(), the vector of returned values
#             has any values below lt replaced by lt, and any values above ut replaced by ut. If
#             resample=yes then no such values should exist anyway.
#
# p.round() This function rounds values to one of the adjacent integers in such a way that no bias
#           is introduced over a large number of samples. For example, if a vector representing
#           10,000 persons has all values set to 4.2, then after p.round is applied, 80% of the
#           people (at random) will get a value of 4, and 20% will get a value of 5.  The overall
#           mean will remain at 4.2.  This is especially important in the product use frequency,
#           for example.  A product may be assigned a point value of (say) 0.4 uses per year, and
#           it would not be useful to round this to zero uses for everybody.  Similarly, a product
#           with a mean usage frequency of 1.4 times per year should get a higher overall usage
#           rate than one with a mean of 0.6 times per year.  Therefore, it is not acceptable to
#           simply round those to once per year and apply that to every person. The p.round function
#           returns two possible values unless the argument is exactly an integer.
#
# quantiles() R contains a function called quantile() already. This one specifies a set of desired
#             quantiles to be evaluated. It is obvious how this set can be modified, if needed.
#
# summarize.chemical() This is called when all the sources have been evaluated for a given chemical.
#                     It is called once per set of persons, and again at the end of the run (if
#                     personOutput=1) to summarize all persons in the run. This version of
#                     summarize.chemical reports statistics on the total population, and on
#                     selected population cohorts, writing all of them to the same file. Each
#                     chemical has its own files.  The statstics are generated by calls to
#                     the summary.stats() function.
#
# summary.stats()  This function is called multiple times for each chemical by summarize.chemical().
#                 summary.stats computes the mean, standard deviation, and various quantiles of the
#                 population distribution for a set of exposure and dose variables, for one chemical
#                 at a time. The same code is used for all population subgroups by subsetting the
#                 output data to the correct group before calling summary.stats().
#
# trimzero()  This function trims extra zeroes off the chemical CAS numbers. THe first part of the
#             CAS number consists of one to six digits, and some databases standardize them by
#             padding with leading zeroes.  Since CAS numbers have to be matched across multiple
#             files in SHEDS, it is necessary to ensure that (for example) 000134_7_4 is matched with
#             134_17_4 on another file, as they would refer to the same chemical (if this was a real
#             CAS number).  If all of the first part consists of zeroes, then trimzero() leaves
#             exactly one, so for example 000000_33_3 becomes 0_33_3.
#
# unpack()  This function is used when ShedsHT is run as an R package. Each time a new working
#           directory is chosen, use unpack() to convert the R data objects into CSV files.
#           Note that both /inputs and /output folders are needed under the chosen directory.
#           Unpack() may be used with a list of object names, in which case just those objects
#           are converted to CSV. This is useful when re-loading one or more defaults into a
#           folder where some of the CSV files have changes, and should not be overwritten.
#           A blank argument or empty list means that all the csv and TXT files in the R
#           package are converted.
#
# vpos()    This function locates an item in a list.  There are probably other R functions that
#           do this, but this one was written because it does not require additional R packages
#           and its behavior can be easily examined.  A faster method might be found, which would
#           make vpos() obZ3solete.
#
# write.persons()  This function writes demographic, exposure, and dose variables to a separate
#                 output file for each chemical.  THis function rounds the variables to a
#                 reasonable precision so that the files are more readable, as unrounded output
#                 is cluttered with entries with (say) 14 digits, most of which are not significant.
#                 Until Nov 2015, the write.persons() functionality was part of the summarize.chemical()
#                 function, but now it is called separately from the main run() function.




