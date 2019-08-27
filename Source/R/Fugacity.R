#' add.fugs
#'
#' Evaluates the variables in the fugacity input file and creates values for those variables corresponding to each simulated person.
#'
#' @param n.per The total number of simulated persons in this model run specified in the \code{\link{Run}} file
#'
#' @param x the output of the \code{\link{read.fug.inputs}} function.
#'
#' @param pdmf he output of the \code{\link{add.factors}} function. A data set containing physiological and demographic parameters for each theoretical person, the duration of exposure to each potential exposure medium for each person,
#'  the media specific exposure factors, and the  number of baths taken and hand wash events occurring per day for each person.
#'
#' @details This function evaluates the variables in the fugacity input file, creating one value for each simulated
#' person and adding each variable as a new column in the pdmf data set (renamed as pdmff). All variables are left i
#' n their original units except for those with units ug/cm2, which are converted to ug/m2. In the fugacity calculations,
#' all masses are in ug and all lengths are in m.
#'
#' @return pdmff Output contains values sampled from the distributions of each relevant variable in the \code{\link{Fugacity}}
#'  input file for each theoretical person.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{add.factors}}, \code{\link{add.media}}, \code{\link{select_people}}, \code{\link{Fugacity}}, \code{\link{run}}
#'
#' @keywords SHEDS
#'
#' @export

add.fugs = function(n.per,x,pdmf) {
  samples <- as.data.frame(matrix(0,nrow=n.per,ncol=nrow(x)))
  for (i in 1:nrow(x)) {
    samples[i] <- distrib(x$form[i],x$par1[i],x$par2[i],x$par3[i],
                          x$par4[i],x$lower.trun[i],x$upper.trun[i],x$resamp[i],n.per)
    setnames(samples,paste0("V",i),x$varname[i])
    if(x$units[i]=="ug/cm2") samples[i] <- samples[i]*1E4
    # 1E4 converts from [ug/cm2] to [ug/m2]
  }
  return(cbind(pdmf,samples))
}

#' chem.fug
#'
#' Creates distributions of chemical specific parameters for each chemical of interest, in order to reflect real-world
#' variability and uncertainty. These distributions are then sampled to create chemical specific parameters associated
#' with the exposure of each simulated person.
#'
#' @param n.per The total number of simulated persons in this model run specified in the \code{\link{Run}} file.
#'
#' @param cprops The chemical properties required for SHEDS-HT (output of \code{chem.props} function). The default data was
#' prepared from publicly available databases using a custom program (not part of SHEDS-HT).
#' The default file contains about 17 numerical inputs per chemical, but most are not used.
#' The required properties are molecular weight (\code{MW}), vapor pressure (\code{VP.Pa}),
#' solubility (\code{water.sol.mg.l}), octanol-water partition coefficient (\code{log.Kow}),
#' air decay rate (\code{half.air.hr}), decay rate on surfaces (\code{half.sediment.hr}),
#' and permeability coefficient (\code{Kp}).
#'
#' @param x From the output of the \code{\link{read.fug.inputs}} function.
#'
#' @details This function obtains chemical-specific properties from the chem.data data set.That data set contains
#' point value for each variable, which this function defines distributions around, with the exception of molecular
#' weight. The constructed distributions reflect both uncertainty and variability (for example, vapor pressure will
#' vary with temperature, humidity, and air pressure or altitude, and may depend on the product formulation).
#' For each variable, a random sample is generated for each simulated person, and the set of variables becomes the
#' output data set.
#' Input surface loading variables are in units of ug/cm2 or ug/cm2/day and are converted by the function to meters,
#' to avoid the need for conversion factors in later equations. These conversions are done after the random sampling,
#' as otherwise the correct conversions depend on the distributional form (for example, par2 is changed for the normal,
#' but not for the lognormal).
#'
#' @return samples A data set with the chemical specific parameters for each combination of chemical and simulated person.
#' For each chemical, the chemical specific parameters assigned to a given person are randomly sampled from distributions
#' on those parameters. These distributions are created from point estimates to reflect real-world uncertainty and variability.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{Fugacity}}, \code{\link{Run}}, \code{\link{run}}, \code{\link{Chem_props}}, \code{\link{get.fug.concs}}
#'
#' @keywords SHEDS
#'
#' @export

chem.fug = function(n.per,cprops,x) {
  x$varname    <- "vapor"
  x$units      <- "Pa"
  x$form       <- "logn"
  if (is.na(cprops$vapor)) { cprops$vapor <- 1
  cat("/n vapor pressure estimated") }
  x$par1       <- min(max(cprops$vapor,1E-12),1E5)
  x$par2       <- 2.0
  x$par3       <- NA
  x$par4       <- NA
  x$lower.trun <- x$par1/100
  x$upper.trun <- min(x$par1*100,1E5)
  x$resamp     <- "yes"
  x$descrip    <- "vapor pressure"
  y            <- x
  x$varname    <- "molwt"
  x$units      <- "g/mol"
  x$form       <- "poin"
  if (is.na(cprops$mw)) { cprops$mw <- 100
  cat("/n molecular weight estimated") }
  x$par1       <- cprops$mw
  x$par2       <- NA
  x$par3       <- NA
  x$par4       <- NA
  x$lower.trun <- NA
  x$upper.trun <- NA
  x$resamp     <- "yes"
  x$descrip    <- "molecular weight"
  y            <- rbind(y,x)
  x$varname    <- "solub"
  x$units      <- "mol/m3"
  x$form       <- "logn"
  if (is.na(cprops$solub)) { cprops$solub <- 1
  cat("/n solubility estimated") }
  x$par1       <- cprops$solub/cprops$mw
  # solub units converted from [mg/l] to [mol/m3]
  # this works before sampling because cprops has point values
  x$par2       <- 1.3
  x$par3       <- NA
  x$par4       <- NA
  x$lower.trun <- x$par1/5
  x$upper.trun <- x$par1*5
  x$resamp     <- "yes"
  x$descrip    <- "solubility"
  y            <- rbind(y,x)
  x$varname    <- "kow"
  x$units      <- "[-]"
  x$form       <- "logn"
  if (is.na(cprops$kow)) { cprops$kow <-10
  cat("/n Kow estimated") }
  x$par1       <- cprops$kow
  x$par2       <- 1.3
  x$par3       <- NA
  x$par4       <- NA
  x$lower.trun <- x$par1/5
  x$upper.trun <- x$par1*5
  x$resamp     <- "yes"
  x$descrip    <- "octanol-water partition coefficient"
  y            <- rbind(y,x)
  x$varname    <- "decay.air"
  x$units      <- "1/day"
  x$form       <- "logn"
  if (is.na(cprops$decay.a)) { cprops$decay.a <- 0.01
  cat("/n air decay rate estimated") }
  x$par1       <- cprops$decay.a
  x$par2       <- 1.5
  x$par3       <- NA
  x$par4       <- NA
  x$lower.trun <- x$par1/10
  x$upper.trun <- x$par1*10
  x$resamp     <- "yes"
  x$descrip    <- "loss rate in air"
  y            <- rbind(y,x)
  x$varname    <- "decay.sur"
  x$units      <- "1/day"
  x$form       <- "logn"
  if (is.na(cprops$decay.s)) { cprops$decay.s <- 0.001
  cat("/n surface decay rate estimated") }
  x$par1       <- cprops$decay.s
  x$par2       <- 1.5
  x$par3       <- NA
  x$par4       <- NA
  x$lower.trun <- x$par1/10
  x$upper.trun <- x$par1*10
  x$resamp     <- "yes"
  x$descrip    <- "loss rate on surfaces"
  y            <- rbind(y,x)
  x$varname    <- "diffus.air"
  x$units      <- "m2/day"
  x$form       <- "logn"
  m            <- cprops$mw
  x$par1       <- 2.05*(1/29 + 1/m)^0.5 / m^0.33
  # from EPA's Chemical Engineering Methods Manual (1991) App. K, units [cm2/s]
  x$par1       <- x$par1 * 86400 / 10000
  # convert from [cm2/s] to [m2/day]
  x$par2       <- 1.3
  x$par3       <- NA
  x$par4       <- NA
  x$lower.trun <- x$par1/5
  x$upper.trun <- x$par1*5
  x$resamp     <- "yes"
  x$descrip    <- "diffusivity in air"
  y            <- rbind(y,x)
  x$varname    <- "h.y0"
  x$units      <- "m/hr"
  x$form       <- "logn"
  x$par1       <- 46.8 * 3.3 / (2.5+cprops$mw^0.333)^2
  x$par2       <- 1.3
  x$par3       <- NA
  x$par4       <- NA
  x$lower.trun <- x$par1/5
  x$upper.trun <- x$par1*5
  x$resamp     <- "yes"
  x$descrip    <- "mass transfer coefficient"
  y            <- rbind(y,x)
  x$varname    <- "c.out.air"
  x$units      <- "ug/m3"
  x$form       <- "poin"
  x$par1       <- 0
  x$par2       <- NA
  x$par3       <- NA
  x$par4       <- NA
  x$lower.trun <- NA
  x$upper.trun <- NA
  x$resamp     <- "yes"
  x$descrip    <- "outdoor chemical concentration in air"
  y <- rbind(y,x)
  x$varname    <- "c.prev.air"
  x$units      <- "ug/m3"
  x$descrip    <- "air concentration before chemical usage"
  y            <- rbind(y,x)
  x$varname    <- "c.prev.sur"
  x$units      <- "ug/cm2"
  x$descrip    <- "surface concentration before chemical usage"
  y            <- rbind(y,x)
  x$varname    <- "c.src.air"
  x$units      <- "ug/m3/day"
  x$descrip    <- "chemical source strength in air"
  y            <- rbind(y,x)
  x$varname    <- "c.src.sur"
  x$units      <- "ug/cm2/day"
  x$descrip    <- "chemical source strength on surfaces"
  y            <- rbind(y,x)
  z            <- as.data.table(y)
  samples <- as.data.frame(matrix(0,nrow=n.per,ncol=nrow(z)))
  for (i in 1:nrow(z)) {
    samples[i] <- distrib(z$form[i],z$par1[i],z$par2[i],z$par3[i],
                          z$par4[i],z$lower.trun[i],z$upper.trun[i],z$resamp[i],n.per)
    setnames(samples,paste0("V",i),z$varname[i])
    if(z$units[i]=="ug/cm2")     samples[i] <- samples[i]*1E4
    if(z$units[i]=="ug/cm2/day") samples[i] <- samples[i]*1E4
    # convert units to [ug/m2] and [ug/m2/day]
  }
  return(samples)
}

#'get.fug.concs
#'
#' Performs fugacity calculations to evaluate time-dependent chemical flows.
#'
#' @param sdata The chemical-scenario data specific to relevant combinations of chemical and scenario. Generated internally.
#'
#' @param chem.data The list of scenario-specific information for the chemicals being evaluated. Generated internally.
#'
#' @param x The output of the \code{\link{add.fugs}} function.
#'
#' @param cfug The output of the \code{\link{chem_fug}} function, a data set with the chemical specific parameters for each
#' combination of chemical and simulated person. For each chemical, the chemical specific parameters assigned to a given person
#' are randomly sampled from distributions on those parameters. These distributions are created from point estimates to reflect
#' real-world uncertainty and variability.
#'
#' @details This is one of two functions that perform fugacity calculations. This one evaluates dynamic or time-dependent chemical
#' First, a set of local variables are determined for use in later calculations.  These are a mix of fixed and chemical-dependent
#' variables (evaluated separately for each person).  Some variables, like chemical mass and app.rates, vary with each source, so
#' these calculations are repeated for each source-scenario.
#' Second, the eigenvalues and eigenvectors of the jacobian matrix are calculated.  Since the fugacity model has been reduced to
#' just two compartments (air and surface), the solutions can be expressed analytically, and there is no explicit invocation of
#' any linear algebra routines that would normally be required.
#' Third, the variables composing the \code{concs} output are evaluated. The variables \code{m.c.air} and \code{m.c.sur}
#' are the time-constant masses, while \code{m.t0.air} and \code{m.t0.sur} are the time-dependent masses at t=0. The
#' time-constant parts are zero here because the permanent sources (i.e. \code{c.src.air} and \code{c.src.sur}) are
#' assumed to be zero in these calculations.  The time-dependent masses are multiplied exponentially as a function of time,
#' and thus approach zero when enough time has passed.
#'
#' @return concs A data set containing calculated dynamic chemical flows for each unique combination
#' of simulated person and chemical.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{chem_fugs}}, \code{\link{Fugacity}}, \code{\link{Run}}, \code{\link{run}}, \code{\link{get.y0.concs}}
#'
#' @keywords SHEDS kwd2
#'
#' @export

get.fug.concs = function(sdata,chem.data,x,cfug) {
  n            <- nrow(x)
  chem.mass    <- 1E6 * sdata$mass * chem.data$f.chemical        # sdata$mass [g], chem.mass [ug]
  app.rate.sur <- chem.mass/x$area.sur                           # app.rate.sur [ug/m2]
  app.rate.air <- app.rate.sur/1E6                               # app.rate.air [ug/m3]
  # above sets init air to same # of ug/m3 as surf conc is in g/m2
  time         <- runif(n) * 365/sdata$use.freq                  # time [days]
  vol.air      <- x$area.sur * x$height                          # vol.air [m3]
  sm.mass.air  <- vol.air * x$sm.load.air                        # sm.mass.air [ug]
  lg.mass.air  <- vol.air * x$lg.load.air                        # lg.mass.air [ug]
  sm.mass.sur  <- x$area.sur * x$sm.load.sur                     # sm.mass.sur [ug]
  lg.mass.sur  <- x$area.sur * x$lg.load.sur                     # lg.mass.sur [ug]
  sm.clean.sur <- pmax(x$sm.clean.sur,x$sm.depos*x$sm.load.air/x$sm.load.sur-x$sm.resus)
  lg.clean.sur <- pmax(x$lg.clean.sur,x$lg.depos*x$lg.load.air/x$lg.load.sur-x$lg.resus)
  # clean.sur [1/day], depos [m/day], load.air [ug/m3], load.sur [ug/m2], resus [1/day]
  sm.clean.air <- pmax(x$sm.clean.air,(x$sm.resus*sm.mass.sur-x$sm.depos*x$sm.load.air*x$area.sur)/sm.mass.air)
  lg.clean.air <- pmax(x$lg.clean.air,(x$lg.resus*lg.mass.sur-x$lg.depos*x$lg.load.air*x$area.sur)/lg.mass.air)
  # clean.air [1/day], resus [1/day], mass.sur [ug], depos [m/day], load.air [ug/m3], area.sur [m2], mass.air [ug]
  ug.mol       <- 1E6*cfug$molwt                                 # ug.mol [ug/mol]
  z.air        <- 1/(8.314*x$temp)                               # z.air [mol/(Pa m3)]
  zvb.air      <- z.air * vol.air * ug.mol                       # zvb.air [ug/Pa]
  z.sur        <- z.air * 82500 / (cfug$vapor^0.65)              # z.sur [mol/(Pa m3)]
  zvb.sur      <- z.sur * x$area.sur * x$thick.sur * ug.mol      # zvb.sur [ug/Pa]
  sm.kp        <- 1.662E-12 * cfug$kow * x$sm.carb.f * cfug$solub / (cfug$vapor * z.air)
  lg.kp        <- 1.662E-12 * cfug$kow * x$lg.carb.f * cfug$solub / (cfug$vapor * z.air)
  # kp [m3/ug], 1.662E-12 [m3/ug], kow [-], carb.f [-], solub [mol/m3], vapor [Pa], z.air [mol/(Pa m3)]
  sm.zv.air    <- zvb.air * sm.kp * x$sm.load.air                # sm.zv.air [ug/Pa]
  lg.zv.air    <- zvb.air * lg.kp * x$lg.load.air                # lg.zv.air [ug/Pa]
  sm.cap       <- z.air * sm.kp * ug.mol                         # sm.cap [1/Pa]
  lg.cap       <- z.air * lg.kp * ug.mol                         # lg.cap [1/Pa]
  sm.zv.sur    <- sm.cap * sm.mass.sur                           # sm.zv.sur [ug/Pa]
  lg.zv.sur    <- lg.cap * lg.mass.sur                           # lg.zv.sur [ug/Pa]
  zv.air       <- zvb.air + sm.zv.air + lg.zv.air                # zv.air [ug/Pa]
  zv.sur       <- zvb.sur + sm.zv.sur + lg.zv.sur                # zv.sur [ug/Pa]
  izv.air      <- pmin(1E100,1/zv.air)                           # izv.air [Pa/ug]
  izv.sur      <- pmin(1E100,1/zv.sur)                           # izv.sur [Pa/ug]
  yaf          <- pmin(cfug$diffus.air*z.air/x$thick.bou , 0.0135/(cfug$vapor^0.32))  # yaf [mol/(m2-Pa-day)]
  cln.air      <- izv.air * (sm.zv.air*sm.clean.air + lg.zv.air*lg.clean.air)         # cln.air [1/day]
  cln.sur      <- izv.sur * (sm.zv.sur*sm.clean.sur + lg.zv.sur*lg.clean.sur)         # cln.sur [1/day]
  sm.dep       <- izv.air * x$area.sur * x$sm.load.air * x$sm.depos * sm.cap          # sm.dep [1/day]
  lg.dep       <- izv.air * x$area.sur * x$lg.load.air * x$lg.depos * lg.cap          # lg.dep [1/day]
  dep          <- sm.dep + lg.dep                                                     # dep [1/day]
  res          <- izv.sur * (sm.mass.sur*x$sm.resus + lg.mass.sur*x$lg.resus)         # res [1/day]
  diff.air     <- izv.air * ug.mol * x$area.sur * yaf                                 # diff.air [1/day]
  diff.sur     <- izv.sur * ug.mol * x$area.sur * yaf                                 # diff.sur [1/day]
  m0.air       <- (cfug$c.prev.air + app.rate.air) * vol.air                          # m0.air [ug]
  m0.sur       <- (cfug$c.prev.sur + app.rate.sur) * x$area.sur                       # m0.sur [ug]
  src.air      <- cfug$c.out.air * x$aer.out * vol.air + cfug$c.src.air * vol.air     # src.air [ug/day]
  src.sur      <- cfug$c.src.sur * x$area.sur                                         # src.sur [ug/day]

  a <- x$aer.out + cfug$decay.air + cln.air + dep + diff.air                          # a [1/day]
  b <- res + diff.sur                                                                 # b [1/day]
  c <- dep + diff.air                                                                 # c [1/day]
  d <- cfug$decay.sur + cln.sur + res + diff.sur                                      # d [1/day]
  r <- sqrt(a^2+4*b*c-2*a*d+d^2)                                                      # r [1/day]
  lam1     <- (a+d+r)/2                                                               # lam1 [1/day]
  lam2     <- (a+d-r)/2                                                               # lam2 [1/day]
  v1.air   <- r+a-d                                                                   # v1.air [1/day]
  v1.sur   <- -2*c                                                                    # v1.sur [1/day]
  v2.air   <- ifelse(r-a+d>0,r-a+d,(2*b*c)/(a-d))                                     # v2.air [1/day]
  v2.sur   <- 2*c                                                                     # v2.sur [1/day]
  m.c.air  <- (d*src.air+b*src.sur)/(a*d-b*c)                                         # m.c.air [ug]
  m.c.sur  <- (c*src.air+a*src.sur)/(a*d-b*c)                                         # m.c.sur [ug]
  m.t0.air <- m0.air - m.c.air                                                        # m.t0.air [ug]
  m.t0.sur <- m0.sur - m.c.sur                                                        # m.t0.sur [ug]
  k1.air   <-  (2*c*m.t0.air+(a-d-r)*m.t0.sur)*v1.air/(4*c*r)                         # k1.air [ug]
  k2.air   <-  (2*c*m.t0.air+(a-d+r)*m.t0.sur)*v2.air/(4*c*r)                         # k2.air [ug]
  k1.sur   <- -(2*c*m.t0.air+(a-d-r)*m.t0.sur)/(2*r)                                  # k1.sur [ug]
  k2.sur   <-  (2*c*m.t0.air+(a-d+r)*m.t0.sur)/(2*r)                                  # k2.sur [ug]

  concs <- data.table(matrix(0,nrow=nrow(cfug),ncol=14))
  setnames(concs,1:14,c("seq","time","mass.air","mass.sur","conc.air","conc.sur",
                        "m0.air","m0.sur","q","k","gain","loss","cq0.air","cq0.sur"))
  concs$seq      <- 1:n                                                               # seq [-]
  concs$time     <- time                                                              # time [day]
  exp1           <- exp(-lam1*time)                                                   # exp1 [-]
  exp2           <- exp(-lam2*time)                                                   # exp2 [-]
  concs$mass.air <- m.c.air + k1.air*exp1 + k2.air*exp2                               # mass.air [ug]
  concs$mass.sur <- m.c.sur + k1.sur*exp1 + k2.sur*exp2                               # mass.sur [ug]
  concs$conc.air <- concs$mass.air / vol.air                                          # conc.air [ug/m3]
  concs$conc.sur <- concs$mass.sur / x$area.sur                                       # conc.sur [ug/m2]
  concs$conc.sur <- ifelse(concs$conc.sur<1E-280,1E-280,concs$conc.sur)
  return(concs)
}

#' get.y0.concs
#'
#' Performs fugacity calculations to evaluate constant (time-independent) chemical flows, such as emissions from household articles.
#'
#' @param sdata The chemical-scenario data specific to relevant combinations of chemical and scenario. Generated internally.
#'
#' @param chem.data The list of scenario-specific information for the chemicals being evaluated. Generated internally.
#'
#' @param pdmff Output from the \code{\link{add.fugs}} function. A data set containing values sampled from the distributions
#' of each relevant variable in the \code{\link{Fugacity}} input file for each theoretical person.
#'
#' @param cfug The output of the \code{\link{chem_fug}} function, a data set with the chemical specific parameters for each
#' combination of chemical and simulated person. For each chemical, the chemical specific parameters assigned to a
#' given person are randomly sampled from distributions on those parameters. These distributions are created from point
#' estimates to reflect real-world uncertainty and variability.
#'
#' @details This function evaluates the chemical concentrations resulting from constant source emissions. Thus, the
#' function employs the steady state solution to the fugacity equations. The basis for these calculations is that the
#' chemical sources are from articles, and are thus permanent and unchanging. The chemical concentrations will therefore
#' quickly adjust so that the flows are balanced, and the concentrations remain fixed thereafter.
#'
#' @return concs A data set containing calculated steady state chemical flows for each unique combination of simulated
#'  person and chemical.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{chem_fugs}}, \code{\link{Fugacity}}, \code{\link{Run}}, \code{\link{run}}, \code{\link{get.y0.concs}}
#'
#' @keywords SHEDS kwd2
#'
#' @export

get.y0.concs = function(sdata,chem.data,pdmff,cfug) {
  z.air       <- 1/(8.314*pdmff$temp)                       # z.air [mol/(Pa m3)]
  sm.kp       <- 1.662E-12 * cfug$kow * pdmff$sm.carb.f * cfug$solub / (cfug$vapor * z.air)
  lg.kp       <- 1.662E-12 * cfug$kow * pdmff$lg.carb.f * cfug$solub / (cfug$vapor * z.air)
  vol.air     <- pdmff$area.sur * pdmff$height             # vol.air [m3]
  # kp [m3/ug], 1.662E-12 [m3/ug], kow [-], carb.f [-], solub [mol/m3], vapor [Pa], z.air [mol/(Pa m3)]
  qstar <- pdmff$aer.out/24 * vol.air * (1 + sm.kp*pdmff$sm.load.air + lg.kp*pdmff$lg.load.air)
  # qstar is Q* in [m3/h], /24 converts aer from [1/day] to [1/hr], vol.air in [m3]
  y           <- chem.data$y0 / (1 + qstar/(cfug$h.y0*sdata$f.area*pdmff$area.sur))
  # y is average gas-phase concentration in house in [ug/m3]
  ug.mol      <- 1E6*cfug$molwt                             # ug.mol [ug/mol]

  zvb.air     <- z.air * vol.air * ug.mol                   # zvb.air [ug/Pa]
  fug         <- y / ug.mol / z.air                         # fug [Pa]
  # fug is the fugacity of the gas-phase component

  conc.air    <- y * ( 1 + sm.kp*pdmff$sm.load.air + lg.kp*pdmff$lg.load.air)
  # conc.air is the (gas+particle) chemical concentration in air [ug/m3]

  mass.air    <- conc.air * vol.air                        # mass.air [ug]
  # mass.air is the airborne chemical mass [ug]
  sm.mass.sur <- pdmff$area.sur * pdmff$sm.load.sur        # sm.mass.sur [ug]
  lg.mass.sur <- pdmff$area.sur * pdmff$lg.load.sur        # lg.mass.sur [ug]
  z.sur       <- z.air * 82500 / (cfug$vapor^0.65)         # z.sur [mol/(Pa m3)]
  zvb.sur     <- z.sur * pdmff$area.sur * pdmff$thick.sur * ug.mol    # zvb.sur [ug/Pa]
  sm.zv.air   <- zvb.air * sm.kp * pdmff$sm.load.air       # sm.zv.air [ug/Pa]
  lg.zv.air   <- zvb.air * lg.kp * pdmff$lg.load.air       # lg.zv.air [ug/Pa]
  sm.cap      <- z.air * sm.kp * ug.mol                    # sm.cap [1/Pa]
  lg.cap      <- z.air * lg.kp * ug.mol                    # lg.cap [1/Pa]
  sm.zv.sur   <- sm.cap * sm.mass.sur                      # sm.zv.sur [ug/Pa]
  lg.zv.sur   <- lg.cap * lg.mass.sur                      # lg.zv.sur [ug/Pa]
  zv.air      <- zvb.air + sm.zv.air + lg.zv.air           # zv.air [ug/Pa]
  zv.sur      <- zvb.sur + sm.zv.sur + lg.zv.sur           # zv.sur [ug/Pa]
  mass.sur    <- fug * zv.sur                              # mass.sur [ug]
  # mass.sur is the chemical mass on or in surfaces
  conc.sur    <- mass.sur / pdmff$area.sur                 # conc.sur [ug/m2]
  # conc.sur is the surface concentration in ug/m2
  # conc.sur    <- ifelse(conc.sur<1E-280,1E-280,conc.sur)
  concs       <- as.data.table(cbind(mass.air,mass.sur,conc.air,conc.sur))
  return(concs)
}

# Notes for SHEDS-HT Fugacity module
#
# add.fugs()  This function evaluates the variables on the fugacity input file,
#             creating one value for each simulated person and adding each variable
#             as a new column on the pdmf data set (renamed as pdmff). All variables
#             are left in their original units except for those with units [ug/cm2],
#             which are converted to [ug/m2]. In the fugacity calculations, all masses
#             are in [ug] and all lengths are in [m].
#
# chen.fug()  This function obtains chemical-specific properities from the chem.data
#             data set, all of which are point values. Here distributions are
#             defined around the reported value for all variables except
#             molecular weight. The distributions reflect both uncertainty and
#             variability (for example, vapor pressure will vary with temperature,
#             humidity, and air pressure or altitude, and may depend on the product
#             formulation).
#
#             Five variables (c.out.air, c.prev.air, c.prev.sur, c.src.air, c.src.sur)
#             are included but set to point values of zero. These could simply be
#             dropped from the resulting calculations, but it would then be more
#             difficult to reinstate them later. If a future version of SHEDS were
#             to allow them to be non-zero, then they would have to be read from an
#             input file, and would have to be chemical-specific (unlike the other
#             variables on the fugacity input file).  For a run with thousands of
#             chemicals, this would become a large input file requiring much effort
#             to research and prepare.
#
#             Surface loading variables using units of [ug/cm2] or [ug/cm2/day] are
#             converted to use standard length units (meters) here, to avoid the need
#             for conversion factors in later equations. These conversions are done
#             after the random sampling, as otherwise the correct conversions depend
#             on the distributional form (for example, par2 is changed for the normal,
#             but not for the lognormal).  For each variable, a random sample is
#             generated for each simulated person, and the set of variables becomes
#             the cfug data set.
#
# get.fug.concs() This is one of two functions that perform fugacity calculations.
#                 This one evaluates dynamic or time-dependent chemical flows. The
#                 function has three parts.
#
#                 First, a set of local variables are determined for use in later
#                 calculations.  These are a mix of fixed and chemical-dependent
#                 variables (evaluated separately for each person).  Some
#                 variables, like chemical mass and app.rates, vary with each
#                 source, so these calculations are repeated for each source-scenario.
#
#                 Second, the eigenvalues and eigenvectors of the jacobian matrix
#                 are calculated.  Since the fugacity model has been reduced to
#                 just two conpartments (air and surface), the solutions can be
#                 expressed analytically, and there is no explicit invocation of
#                 any linear algebra routines that would normally be required.
#
#                 Third, the variables in the output "concs" data set are evaluated.
#                 The variables m.c.air and m.c.sur are the time-constant masses,
#                 while m.t0.air and m.t0.sur are the time-dependent masses at t=0.
#                 The time-constant parts are zero in the absence of permanent sources
#                 (i.e. c.src.air and c.src.sur), so they should always be zero here
#                 because those sources are set to zero.  The time-dependent masses
#                 are multiplied by exponentials in time to obtain chemical masses
#                 at later times, so eventually these approach zero when enough time
#                 has passed.
#
#
# get.y0.concs()  This function evaluates the chemical concentrations in the y0
#                 scenario. This uses the steady state solution to the fugacity
#                 equations. The argument is that the chemical sources are from
#                 articles, and are permanent and unchanging. The chemical
#                 concentrations will quickly adjust so that the flows are
#                 balanced, and the concentrations remain fixed thereafter.
#                 While it is true that the chemical emission rates may decrease
#                 as an article ages, at some point the azricle will wear out and
#                 be replaced with a new one.  SHEDS-HT is intended as a low tier
#                 screening model, and time-varying article concentrations are not
#                 considered necessary at this time.


