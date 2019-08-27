# This is the main module of EPA's SHEDS-HT chemical exposure model.
# This model is derived from SHEDS-Lite; before that SHEDS-Multimedia.
# First R version programmed by WGG in 2012
# Latest changes by KKI on Nov 16, 2017
# To run SHEDS-HT, open this module in R or R-studio, source this module
# and type run() to use the default run file "Run_test.txt",
# or specify an alternate run file (in quotes) as an argument to run().
# All input files must be in the /inputs folder under the working
# directory wd. If moving this code to a new location, change the
# default "wd=" argument in the run() function call immediately below.

#' run
#'
#' Function to call the \code{\link{Run}} txt file, which consists of user-defined parameters and calls to input files required
#' to initialize a SHEDS.HT run.
#'
#' @param run.file The name of the run file to be used for a given run.  Many different "Run"" files may be set up for special
#' purposes.The one being invoked must be present in the inputs folder.
#'
#' @param wd The user's working directory. The working directory should contain an inputs folder, containing all necessary
#' SHEDS.HT input files. The wd value should be set once for each SHEDS.HT installation by replacing the default value in
#' the definition of run().
#'
#' @details This function is used to provide R with information necessary to produce a SHEDS-HT run and to initialize necessary
#' parameters. No values are returned. In order to produce a successful run, all necessary input files should be stored in the
#'  working directory specified by the wd argument.
#'
#' @return No variable will be returned.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{Run.txt}}, \code{\link{read.run.file}}
#'
#' @keyword kwd1
#'
#' @export

run = function(run.file="", wd="") {
  setup(wd)
  specs         <- read.run.file(run.file)
  act.diaries   <- read.act.diaries(specs$act.diary.file,specs)
  chem.props    <- read.chem.props(specs$chem.props.file,specs)
  specs         <- update.specs(specs,chem.props)
  diet.diaries  <- read.diet.diaries(specs$diet.diary.file,specs)
  exp.factors   <- read.exp.factors(specs$exp.factor.file)
  fug.vars      <- read.fug.inputs(specs$fugacity.file)
  media         <- read.media.file(specs$media.file)
  media.sur     <- tolower(media$media[media$type!="air"])
  media.air     <- tolower(media$media[media$type=="air"])
  physio        <- read.phys.file(specs$physiology.file)
  pop           <- read.pop.file(specs$population.file,specs)
  source.scen   <- read.source.scen.file(specs$source.scen.file)
  source.chem   <- read.source.chem.file(specs$source.chem.file, source.scen$src,specs)
  specs         <- update.specs(specs,source.chem)
  source.scen   <- source.scen[source.scen$src %in% source.chem$src]
  source.vars   <- read.source.vars.file(specs$source.vars.file,source.scen)

  act.pools     <- act.diary.pools(act.diaries,specs)
  diet.pools    <- diet.diary.pools(diet.diaries,specs)
  gen.facs      <- gen.factor.tables(exp.factors)
  med.facs      <- med.factor.tables(exp.factors,media.sur)
  sv            <- set.pars(source.vars)
  scv           <- set.pars(source.chem)

  sets <- ceiling(specs$n.persons/specs$set.size)

  for (set in 1:sets) {
    if (set<sets) n.per <- specs$set.size else
      n.per <- specs$n.persons-specs$set.size*(set-1)
    pd    <- select.people(n.per,pop,physio,act.pools,diet.pools,
                           act.diaries,diet.diaries,specs)
    pdm   <- add.media(n.per,media,pd)
    pdmf  <- add.factors(n.per,gen.facs,med.facs,exp.factors,media.sur,pdm)
    pdmff <- add.fugs(n.per,fug.vars,pdmf)
    base  <- pdmff[order(person)]
    src.list  <- source.scen$src
    n.src     <- length(src.list)
    svar.list <- unique(source.vars$varname)
    n.svar    <- length(svar.list)
    src.data  <- array(0,c(n.src,n.per,n.svar))
    n.sdist   <- nrow(sv)
    if (n.sdist>0) {
      for (s in 1:n.sdist) {
        q <- runif(n.per)
        x <- distrib(sv[s]$form,sv[s]$par1,sv[s]$par2,sv[s]$par3,sv[s]$par4,
                     sv[s]$lower.trun,sv[s]$upper.trun,'y',n=n.per,q=q,v=sv[s]$values)
        y <- pd$age>=sv[s]$min.age & pd$age<=sv[s]$max.age
        if(sv[s]$gender!="") y <- y & sv[s]$gender==pd$gender
        src.data[vpos(sv[s]$src,src.list),y,vpos(sv[s]$varname,svar.list)] <- x[y]
      } }
    for (c in 1:specs$n.chem) {
      if(exists("fexp")) rm(fexp, inherits = TRUE)

      chem      <- specs$chem.list[c]
      cb        <- make.cbase(base,chem)
      cprops    <- chem.props[chem.props$cas==chem]
      cfug      <- chem.fug(n.per,cprops,fug.vars[1])
      schem     <- scv[scv$cas==chem]
      sources   <- unique(schem$src)
      n.csrc    <- length(sources)
      cvar.list <- unique(schem$varname)
      n.cvar    <- length(cvar.list)
      chem.data <- as.data.frame(array(0,c(n.per,n.cvar)))
      setnames(chem.data,names(chem.data),cvar.list)
      src.means <- as.data.table(array(0,c(n.csrc+1,8)))
      setnames(src.means,names(src.means),c("exp.dermal","exp.ingest","exp.inhal","dose.inhal",
                                            "f.dermal","f.ingest","f.inhal","mean.mass"))
      src.names <- c(sources,"Total")
      for (s in 1:n.csrc) {
        cat("\n Starting source ",s, " of chem ",c,"(",cprops$cas,") of ",specs$n.chem)
        src   <- sources[s]
        scsel <- schem$src==src
        sc    <- schem[scsel]
        for (v in 1:nrow(sc)) {
          q <- runif(n.per)
          chem.data[vpos(sc[v]$varname,cvar.list)] <- distrib(sc[v]$form, sc[v]$par1,
                                                              sc[v]$par2, sc[v]$par3, sc[v]$par4, sc[v]$lower.trun,
                                                              sc[v]$upper.trun,'y',n=n.per,q=q,v=sc[v]$values)
        }
        if(exists("f.chemical",chem.data)) mode(chem.data$f.chemical)<-"numeric"
        sel   <- source.scen$src==src
        scens <- source.scen[sel]
        sdata <- as.data.table(src.data[vpos(src,src.list),,])
        if (length(svar.list>0)) setnames(sdata,names(sdata),svar.list)
        io <- scens$indoor
        if (any(names(sdata)=="use.freq")) {sdata$use.today<-p.round(sdata$use.freq/365)}
        if(scens$dietary==1) {
          dietary <- food.residue(chem.data,cb,src)
          add     <- dietary$exp.ingest.dietary
          src.means$exp.ingest[s] <- src.means$exp.ingest[s] + mean(add)
          src.means$f.ingest[s]   <- 1
          cb$exp.ingest.tot  <- cb$exp.ingest.tot  + add
          cb$exp.dietary.tot <- cb$exp.dietary.tot + add
        }
        if(scens$dirderm==1) {
          dir.derm                <- dir.dermal(sdata,chem.data)
          add                     <- dir.derm$exp.dermal.dirderm
          src.means$exp.dermal[s] <- src.means$exp.dermal[s] +mean(add)
          mean.mass               <- mean(sdata$mass*chem.data$f.chemical*chem.data$chem.prev*1E6)
          src.means$mean.mass[s]  <- mean.mass
          if(mean.mass>0) src.means$f.dermal[s] <- src.means$exp.dermal[s]/mean.mass
          cb$exp.dermal.tot       <- cb$exp.dermal.tot + add
        }
        if(scens$diringest==1) {
          dir.ingest              <- dir.ingested(sdata,chem.data)
          add                     <- dir.ingest$exp.ingest.diringest
          src.means$exp.ingest[s] <- src.means$exp.ingest[s] + mean(add)
          mean.mass               <- mean(sdata$mass*chem.data$f.chemical*chem.data$chem.prev*1E6)
          src.means$mean.mass[s]  <- mean.mass
          if(mean.mass>0) src.means$f.ingest[s] <- src.means$exp.ingest[s]/mean.mass
          cb$exp.ingest.tot       <- cb$exp.ingest.tot  + add
          cb$exp.nondiet.tot      <- cb$exp.nondiet.tot + add
        }
        if(scens$dirinhaer==1) {
          dir.inh.aer             <- dir.inhal.aer(sdata,chem.data,cb,io)
          addexp                  <- dir.inh.aer$exp.inhal.dirinhaer
          adddose                 <- dir.inh.aer$dose.inhal.dirinhaer
          src.means$exp.inhal[s]  <- src.means$exp.inhal[s]  + mean(addexp)
          src.means$dose.inhal[s] <- src.means$dose.inhal[s] + mean(adddose)
          mean.mass               <- mean(sdata$mass*chem.data$f.chemical*chem.data$chem.prev*1E6)
          src.means$mean.mass[s]  <- mean.mass
          if(mean.mass>0) src.means$f.inhal[s] <- src.means$dose.inhal[s]/mean.mass
          cb$exp.inhal.tot        <- cb$exp.inhal.tot  + addexp
          cb$conc.inhal.max.prod.aer  <- pmax(cb$conc.inhal.max.prod.aer, dir.inh.aer$conc)
          cb$dose.inhal.tot       <- cb$dose.inhal.tot + adddose
        }
        if(scens$dirinhvap==1) {
          dir.inh.vap             <- dir.inhal.vap(sdata,chem.data,cprops,cb,io)
          addexp                  <- dir.inh.vap$exp.inhal.dirinhvap
          adddose                 <- dir.inh.vap$dose.inhal.dirinhvap
          src.means$exp.inhal[s]  <- src.means$exp.inhal[s]  + mean(addexp)
          src.means$dose.inhal[s] <- src.means$dose.inhal[s] + mean(adddose)
          mean.mass               <- mean(sdata$mass*chem.data$f.chemical*chem.data$chem.prev*1E6)
          src.means$mean.mass[s]  <- mean.mass
          if(mean.mass>0) src.means$f.inhal[s] <- src.means$dose.inhal[s]/mean.mass
          cb$exp.inhal.tot  <- cb$exp.inhal.tot  + addexp
          cb$conc.inhal.max.prod.vap  <- pmax(cb$conc.inhal.max.prod.vap, dir.inh.vap$conc)
          cb$dose.inhal.tot <- cb$dose.inhal.tot + adddose
        }
        if(scens$downthedrain==1) {
          exp.ddd.mass   <- down.the.drain.mass(sdata,chem.data)
          cb$exp.ddd.tot <- cb$exp.ddd.tot + exp.ddd.mass
        }
        if(scens$indir.fug==1) {
          fug.concs   <- get.fug.concs(sdata,chem.data,pdmff,cfug)
          outwindow   <- fug.concs$mass.air * pdmff$aer.out
          indir       <- indir.exposure(sdata,cb,fug.concs,chem.data)
          src.means$exp.dermal[s] <- src.means$exp.dermal[s] + mean(indir$exp.dermal.indirect)
          src.means$exp.ingest[s] <- src.means$exp.ingest[s] + mean(indir$exp.ingest.indirect)
          src.means$exp.inhal[s]  <- src.means$exp.inhal[s]  + mean(indir$exp.inhal.indirect)
          src.means$dose.inhal[s] <- src.means$dose.inhal[s] + mean(indir$dose.inhal.indirect)
          mean.mass               <- mean(sdata$mass*chem.data$chem.prev*1E6)
          src.means$mean.mass[s]  <- mean.mass
          if(mean.mass>0) {
            src.means$f.dermal[s] <- src.means$exp.dermal[s]/mean.mass
            src.means$f.ingest[s] <- src.means$exp.ingest[s]/mean.mass
            src.means$f.inhal[s]  <- src.means$dose.inhal[s]/mean.mass
          }
          cb$exp.dermal.tot  <- cb$exp.dermal.tot  + indir$exp.dermal.indirect
          cb$exp.ingest.tot  <- cb$exp.ingest.tot  + indir$exp.ingest.indirect
          cb$exp.nondiet.tot <- cb$exp.nondiet.tot + indir$exp.ingest.indirect
          cb$exp.inhal.indir        <- cb$exp.inhal.indir  + indir$exp.inhal.indirect
          cb$exp.inhal.tot   <- cb$exp.inhal.tot   + indir$exp.inhal.indirect
          cb$dose.inhal.tot  <- cb$dose.inhal.tot  + indir$dose.inhal.indirect
          cb$exp.window.tot  <- cb$exp.window.tot  + outwindow
        }
        if(scens$indir.y0==1) {
          y0.concs    <- get.y0.concs(sdata,chem.data,pdmff,cfug)
          outwindow   <- y0.concs$mass.air * pdmff$aer.out
          indir       <- indir.exposure(sdata,cb,y0.concs,chem.data)
          src.means$exp.dermal[s] <- src.means$exp.dermal[s] + mean(indir$exp.dermal.indirect)
          src.means$exp.ingest[s] <- src.means$exp.ingest[s] + mean(indir$exp.ingest.indirect)
          src.means$exp.inhal[s]  <- src.means$exp.inhal[s]  + mean(indir$exp.inhal.indirect)
          src.means$dose.inhal[s] <- src.means$dose.inhal[s] + mean(indir$dose.inhal.indirect)
          y0.tot.mass             <- y0.concs$mass.air + y0.concs$mass.sur
          mean.mass               <- mean(y0.tot.mass*chem.data$chem.prev*1E6)
          src.means$mean.mass[s]  <- mean.mass
          if(mean.mass>0) {
            src.means$f.dermal[s] <- src.means$exp.dermal[s]/mean.mass
            src.means$f.ingest[s] <- src.means$exp.ingest[s]/mean.mass
            src.means$f.inhal[s]  <- src.means$dose.inhal[s]/mean.mass
          }
          cb$exp.dermal.tot  <- cb$exp.dermal.tot  + indir$exp.dermal.indirect
          cb$exp.ingest.tot  <- cb$exp.ingest.tot  + indir$exp.ingest.indirect
          cb$exp.nondiet.tot <- cb$exp.nondiet.tot + indir$exp.ingest.indirect
          cb$exp.inhal.tot   <- cb$exp.inhal.tot   + indir$exp.inhal.indirect
          cb$dose.inhal.tot  <- cb$dose.inhal.tot  + indir$dose.inhal.indirect
          cb$exp.window.tot  <- cb$exp.window.tot  + outwindow
        }
        if(scens$migration==1){
          migration <- food.migration(chem.data,sdata,cb,src)
          add       <- migration$exp.ingest.migrat
          src.means$exp.ingest[s] <- src.means$exp.ingest[s] + mean(add)
          src.means$f.ingest[s]   <- 1
          cb$exp.ingest.tot <- cb$exp.ingest.tot  + add
          cb$exp.migrat.tot <- cb$exp.migrat.tot + add
        }
      } # end loop over sources
      fexp   <- post.exposure(cb,cprops)
      summarize.chemical(fexp,c,chem,cprops$chem.name,set,sets,specs)
      if (specs$person.output==1) write.persons(fexp,chem,set,specs)
      if (specs$source.output==1) {
        dir  <- paste0("output/",specs$run.name)
        name <- paste0(dir,"/CAS_",chem,"_set",set,"_srcMeans.csv")
        src.means[n.csrc+1]$exp.dermal <- sum(src.means$exp.dermal)
        src.means[n.csrc+1]$exp.ingest <- sum(src.means$exp.ingest)
        src.means[n.csrc+1]$exp.inhal  <- sum(src.means$exp.inhal)
        src.means[n.csrc+1]$dose.inhal <- sum(src.means$dose.inhal)
        write.csv(cbind(src.names,src.means),name)
      }
    } # end loop over chemical
  } # end loop over sets

  if (specs$person.output==1) {
    dir  <- paste0("output/",specs$run.name)
    for (c in 1:specs$n.chem) {
      chem    <- specs$chem.list[c]
      name    <- paste0(dir,"/CAS_",chem,"_all.csv")
      cprops  <- chem.props[chem.props$cas==chem]
      results <- as.data.table(fread(name))
      summarize.chemical(results,c,chem,cprops$chem.name,"allstats",sets,specs)
    }
  }
  if (specs$source.output==1) {
    direc  <- paste0("output/",specs$run.name)
    for (c in 1:specs$n.chem) {
      chem    <- specs$chem.list[c]
      for (set in 1:sets) {
        name <- paste0(direc,"/CAS_",chem,"_set",set,"_srcMeans.csv")
        tmp  <- read.csv(name)
        tmp[1] <- NULL
        wgt  <- specs$set.size/specs$n.persons
        if (set==sets) wgt <- 1 - specs$set.size*(sets-1)/specs$n.persons
        if (set==1) { accum <- data.frame(array(0,c(nrow(tmp),9)))
        names(accum)   <- names(tmp)
        accum$src.names <- tmp$src.names
        }
        accum$exp.dermal <- accum$exp.dermal + tmp$exp.dermal * wgt
        accum$exp.ingest <- accum$exp.ingest + tmp$exp.ingest * wgt
        accum$exp.inhal  <- accum$exp.inhal  + tmp$exp.inhal  * wgt
        accum$dose.inhal <- accum$dose.inhal + tmp$dose.inhal * wgt
        accum$f.dermal   <- accum$f.dermal   + tmp$f.dermal   * wgt
        accum$f.ingest   <- accum$f.ingest   + tmp$f.ingest   * wgt
        accum$f.inhal    <- accum$f.inhal    + tmp$f.inhal    * wgt
        accum$mean.mass  <- accum$mean.mass  + tmp$mean.mass  * wgt
      }  # end loop over sets
      tot.mean.mass <- sum(accum$mean.mass)
      tot.f.dermal  <- sum(accum$f.dermal*accum$mean.mass)
      tot.f.ingest  <- sum(accum$f.ingest*accum$mean.mass)
      tot.f.inhal   <- sum(accum$f.inhal *accum$mean.mass)
      x <- as.data.table(accum)
      y <- nrow(x)
      x$mean.mass[y] <- tot.mean.mass
      x$f.dermal[y]  <- tot.f.dermal/tot.mean.mass
      x$f.ingest[y]  <- tot.f.ingest/tot.mean.mass
      x$f.inhal[y]   <- tot.f.inhal/ tot.mean.mass
      name <- paste0(dir,"/CAS_",chem,"_all_srcMeans.csv")
      write.csv(x,name)
    } # end loop over chemicals
  } # end if block on source.output
} # end run

#' setup
#'
#' Loads required R packages and sources the modules needed to perform a SHEDS.HT run. The user might need to
#' download the packages if they are not already present (see Dependencies).
#'
#' @param wd The User's working directory (set in the inputs to the \code{\link{run}} function). Should contain
#' all necessary SHEDS.HT inputs in an /Input folder.
#'
#' @return No values returned.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @note Requires: \code{\link{data.table}}, \code{\link{plyr}}, \code{\link{stringr}}, \code{\link{ggplot2}}
#'
#' @seealso \code{\link{run}}
#'
#' @keyword SHEDS.HT
#'
#' @export

setup = function(wd="") {
  # WD is the working directory, assumed to be one level up from the R,
  # input, and output directories. File path names are relative to this.
  #Disclaimer
  cat("\nShedsHT Version 0.1.7 (02/21/2019)")
  cat("\nDisclaimer")
  cat("\nThe United States Environmental Protection Agency through its Office of Research and Development
  funded and collaborated in the research and development of this software, in part under Contract EP-C-14-001
  to ICF International. The model is publicly available in Beta version form. All input data used for a given
  application should be reviewed by the researcher so that the model results are based on appropriate data
  sources for the given application. This model, default input files, and R package are under continued development
  and testing. The model equations and approach are published in the peer-reviewed literature
  (Isaacs et al. Environ. Sci. Technol. 2014, 48, 12750-12759). The data included herein do not represent
  and should not be construed to represent any Agency determination or policy.")
  if (wd!="") setwd(wd)
  suppressPackageStartupMessages(TRUE)
  # Load required packages
  library("data.table")
  library("stringr")
  library("plyr")
  library("ggplot2")
}

#' act.diary.pools
#'
#' Assigns activity diaries from the \code{\link{Activity_diaries}} input on the \code{\link{Run}} file
#' (read through the \code{\link{read.act.diaries}} function) to pools based on age, gender, and season
#'
#' @param diaries A data set created internally in SHEDS.HT through the \code{\link{read.act.diaries}} function. The data are
#' activity diaries, which indicate the amount of time and level of metabolic activity in various micro environments.
#' Each line of data represents one person-day (24 hours).
#'
#' @param specs Output of the \code{\link{read.run.file}} function, which can be modified by the \code{\link{update.specs}}
#' function before input into \code{act.diary.pools}.
#'
#' @details The \code{act.diary.pools} function assigns the activity diaries to pools. A given diary may belong to
#' many pools, since every year of age has its own pool. Large pools may contain a list of several hundred diary numbers.
#' The code contains four loops and on each step performs a sub-setting of the list of diary numbers.
#'
#' @return pool A vector of lists.  The length of \code{pool} is the product of the genders, seasons, and ages inputs specified
#' in the \code{\link{Run}} file. Each element is a list of acceptable activity diary numbers for each year of age, gender,
#' weekend, and season combination. For example, \code{pool}[100] may have the name M0P99, which indicates that it is for
#' males, on weekdays, in spring,for age=99.
#' In addition, If the function runs successfully, the following message will be printed:
#' "Activity Diary Pooling completed"
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @note The input to \code{act.diary.pools}, diaries, is created by reading in the Activity_diaries input file (specified on
#' the \code{\link{Run}} file) with the \code{\link{read.act.diaries}} function within the \code{\link{run}} function.
#'
#' @seealso \code{\link{run}}
#'
#' @keyword SHEDS
#'
#' @export

act.diary.pools = function(diaries,specs) {
  i <- 0
  pool <- vector("list",length(specs$genders)*2*
                   length(specs$seasons)*(1+specs$max.age-specs$min.age))
  for (g in specs$genders) {
    for (w in c(0,1)) {
      for (s in specs$seasons) {
        for (age in specs$min.age:specs$max.age) {
          i <- i+1
          agediff <- round(age*specs$age.match.pct/100)
          if (specs$age.match.pct>0) agediff <- max(agediff,1)
          names(pool)[[i]] <- paste0(g,w,s,age)
          pool[i] <- list(diaries$diary.id[diaries$d.gender==g &
                                             diaries$d.weekend==w & diaries$d.season==s &
                                             diaries$d.age>=age-agediff & diaries$d.age<=age+agediff])
        }
      }
    }
  }
  cat("\n Activity Diary Pooling completed" )
  return(pool)
}

#' diet.diary.pools
#'
#' Assigns activity diaries from the \code{\link{Diet_diaries}} input on the \code{\link{Run}} file
#' (read through the \code{\link{read.run.file}} function) to pools based on age and gender.
#'
#' @param diet.diaries A data set created internally in SHEDS.HT through the \code{\link{read.diet.diaries}} function.
#' The data are daily diaries of dietary consumption by food group. Each line represents one person-day, with demographic
#' variables followed by amounts (in grams/day) for a list of food types indicated by a short abbreviation on the header line.
#'
#' @param specs Output of the \code{\link{read.run.file}} function, which can be modified by the \code{\link{update.specs}}
#' function before input into \code{diet.diary.pools}.
#'
#' @details The \code{diet.diary.pools} function assigns the dietary diaries to pools. A given diary may belong to many pools,
#' since every year of age has its own pool. Large pools may contain a list of several hundred diary numbers. The code contains
#' four loops and on each step performs a sub-setting of the list of diary numbers.
#'
#' @return dpool A vector of lists.  The length of \code{dpool} is the product of the gender  and age inputs specified in
#' the \code{\link{Run}} file. Each element is a list of acceptable diet diary numbers for each year of age and each gender
#' combination.
#' In addition, if the function runs successfully, the following message will be printed: "Dietary Diary Pooling completed"
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @note The input to \code{diet.diary.pools}, \code{diet.diaries}, is created by reading in the \code{\link{Diet_diaries}}
#' input file (specified on the \code{\link{Run}} file) with the \code{\link{read.diet.diaries}} function within the
#' \code{\link{run}} function.
#'
#' @seealso \code{\link{Diet_diaries}}, \code{\link{run}}, \code{\link{read.run.file}}, \code{\link{Run}}
#'
#' @keyword SHEDS
#'
#' @export

diet.diary.pools = function(diaries,specs) {
  i <- 0
  dpool <- vector("list",length(specs$genders)*(1+specs$max.age-specs$min.age))
  for (g in specs$genders) {
    for (age in specs$min.age:specs$max.age) {
      i <- i+1
      agediff <- round(age*specs$age.match.pct/100)
      if (specs$age.match.pct>0) agediff <- max(agediff,1)
      names(dpool)[[i]] <- paste0(g,age)
      dpool[i] <- list(diaries$diet.id[diaries$f.gender==g &
                                         diaries$f.age>=age-agediff & diaries$f.age<=age+agediff])
    }
  }
  cat("\n Dietary Diary Pooling completed")
  return(dpool)
}

#' gen.factor.tables
#'
#' Constructs tables of non-media specific exposure factors for each relevant combination of age, gender, and season. The
#' input is from the \code{\link{Exp_actors}} file (specified on the \code{\link{Run}} file) after being read through the
#' \code{\link{read.exp.factors}} function.
#'
#' @param ef A data set created internally using the \code{\link{run}} function and the \code{\link{read.exp.factors}}
#' function to import the user specified \code{\link{Exp_factors}} input file. The data set contains the distributional
#' parameters for the exposure factors. All of these variables may have age or gender-dependent distributions, although
#' in the absence of data, many are assigned a single distribution from which all persons are sampled.
#'
#' @return exp.gen Output consists of non-media specific exposure factors as a data table. Depending on the user input,
#' these generated exposure factors may be gender and/or season specific.  Each gender-season combination is assigned
#' a row number pointing to the appropriate distribution on the output dataset from the \code{\link{read.exp.factors}} function.
#' Thus, \code{exp.gen} consists of 8 rows per variable (2 genders x 4 seasons), regardless of the number of different
#' distributions used for that variable.
#' In addition, if the function runs successfully, the following message is printed: "General Factor Tables completed"
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @note The input to \code{gen.factor.tables} is created by reading in the \code{\link{Exp_factors}} input file specified on
#' the \code{\link{Run} file with the \code{\link{read.exp.factors}} function within the \code{\link{run}} function.
#'
#' @seealso \code{\link{Exp_factors}}, \code{\link{Run}}, \code{\link{run}}, \code{\link{read.exp.factors}}
#'
#' @keyword SHEDS
#'
#' @export

gen.factor.tables = function(ef=exp.factors) {
  media.facs <- c("avail.f","dermal.tc","om.ratio")
  vars<-c("varname","min.age","max.age","gender","season","media","row")
  dt <- as.data.table(as.data.frame(ef)[vars])
  expgen <- dt
  rgen <- nrow(dt)
  for (i in 1:nrow(dt)) {
    if (!dt$varname[i] %in% media.facs) {
      if (is.na(dt$gender[i])|dt$gender[i]=="") { g <- c("F","M")
      } else g <- as.character(dt$gender[i])
      if (is.na(dt$season[i])|dt$season[i]=="") { s <- c("W","P","S","F")
      } else s <- as.character(dt$season[i])
      for (ig in 1:length(g)) {
        gend <- g[[ig]]
        for (is in 1:length(s)) {
          seas <- s[[is]]
          rgen <- rgen+1
          expgen <- rbind(expgen,dt[i,])
          expgen$gender[[rgen]] <- gend
          expgen$season[[rgen]] <- seas
        }
      }
    }
  }
  exp.gen <- expgen[(nrow(dt)+1):nrow(expgen)]
  cat("\n General Factor Tables completed")
  return(exp.gen)
}

#' med.factor.tables
#'
#' Constructs tables of media specific exposure factors for each relevant combination of age, gender, and season, and each
#' of three media specific variables in the \code{ef} argument. These are avail.f, dermal.tc, and om.ratio.
#'
#' @param ef A data set created internally using the \code{\link{run}} function and the \code{\link{read.exp.factors}} function
#' to import the user specified \code{\link{Exp_factors}} input file. The data set contains the distributional parameters for the
#' exposure factors. All of these variables may have age or gender-dependent distributions, although in the absence of data, many
#' are assigned a single distribution from which all persons are sampled.
#'
#' @param media.sur A list of surface media. This data set is created internally by sub-setting the \code{\link{Media}} input
#' file (read in with the \code{\link{read.media.file}} function within the \code{\link{run}} function) to extract only surface
#' media.
#'
#' @details The three media specific variables in the \code{ef} argument are as follows:
#' \code{avail.f}{Fraction of chemical available for transfer from surfaces via touching.}
#' \code{dermal.tc}}{Dermal transfer coefficient, cm^2/hr.}
#' \code{om.ratio}}{Ratio of object-to-mouth exposure to indirect dermal exposure.}
#'
#' @return exp.med Media specific exposure factors presented as a data table. The present version of the model consists of only
#' 3 such exposure factors, but the code will accept more. Depending on the user input, these generated exposure factors may be
#' gender and/or season specific in addition to media specific. The output consists of 24 rows per variable (2 genders x 4
#' seasons x 3 media), when all ages share the same distribution. If a variable has N age categories (each with its own
#' distribution) then there are (24 x N) rows for that variable.
#' In addition, if the function runs successfully, the following message is printed: "Media-specific Factor Tables completed"
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @note The first input argument to \code{med.factor.tables} (ef) is created by reading in the \code{\link{Exp_factors}}
#' input file (specified on the \code{\link{Run}} file) with the \code{\link{read.exp.factors}} function within the
#' \code{\link{run}} function. The \code{media.sur} input is created by sub-setting the \code{\link{Media}} input file (read in
#' with the \code{\link{read.media.file}} function within the \code{\link{run}} function).
#'
#' @seealso \code{\link{Media}}, \code{\link{Exp_factors}}, \code{\link{Run}}, \code{\link{run}}, \code{\link{read.exp.factors}}, \code{\link{read.media.file}}
#'
#' @keyword SHEDS
#'
#' @export

med.factor.tables = function(ef, media.sur) {
  media.facs <- c("avail.f","dermal.tc","om.ratio")
  vars<-c("varname","min.age","max.age","gender","season","media","row")
  df <- as.data.frame(ef)[vars]
  dt <- as.data.table(df)
  expmed <- dt
  rmed <- nrow(dt)
  for (i in 1:nrow(dt)) {
    if (dt$varname[i] %in% media.facs) {
      if (is.na(dt$gender[i])|dt$gender[i]=="") { g <- c("F","M")
      } else g <- as.character(dt$gender[i])
      if (is.na(dt$season[i])|dt$season[i]=="") { s <- c("W","P","S","F")
      } else s <- as.character(dt$season[i])
      if (is.na(dt$media[i])|dt$media[i]=="") { m <- media.sur
      } else m <- tolower(as.character(dt$media[i]))
      for (ig in 1:length(g)) {
        gend <- g[[ig]]
        for (is in 1:length(s)) {
          seas <- s[[is]]
          for (im in 1:length(m)) {
            med <- m[[im]]
            rmed <- rmed+1
            expmed <- rbind(expmed,dt[i,])
            expmed$gender[[rmed]] <- gend
            expmed$season[[rmed]] <- seas
            expmed$media[[rmed]]  <- med
          }
        }
      }
    }
  }
  exp.med <- expmed[(nrow(dt)+1):nrow(expmed)]
  cat("\n Media-specific Factor Tables completed \n")
  return(exp.med)
}

#' set.pars
#' Adjusts certain SHEDS input distribution types to conform with the requirements of the \code{\link{distrib}} function.
#'
#' @param vars Output of the \code{\link{read.source.chem.file}} or the \code{\link{read.source.vars.file}} functions.
#'
#' @details Set.pars adjusts certain SHEDS input distribution types to conform with the requirements of the Distrib function.
#' The lognormal parameters are changed from arithmetic mean and standard deviation to geometric mean and geometric standard
#' deviation.  For the normal distribution, the standard deviation is computed from the mean and coefficient of variation (CV).
#' For user prevalence, the input may be specified either as a point value (indicating probability) or as a Bernoulli distribution.
#' If the former is used, it is converted to the latter.
#'
#' @return v The modified version of the input vars.
#'
#'@author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{run}}, \code{\link{read.source.chem.file}}, \code{\link{read.source.vars.file}}
#'
#' @keyword SHEDS
#'
#' @export

set.pars = function(vars) {
  v     <- vars
  v$row <- 1:nrow(v)
  v$form[v$varname=="use.prev"  & v$form=="point"] <- "bernoulli"
  v$form[v$varname=="home.prev" & v$form=="point"] <- "bernoulli"
  v$form[v$form=="lognormal" & v$mean==0] <- "point"
  f <- strtrim(tolower(v$form),4)
  v$par1 <- v$mean
  v$par1[f=="logn"] <- v$mean[f=="logn"] /(1+v$cv[f=="logn"]^2)^0.5
  v$par2 <- NA
  v$par2[f=="logn"] <- exp(log(1+v$cv[f=="logn"]^2)^0.5)
  v$par2[f=="norm"] <- v$mean[f=="norm"]*v$cv[f=="norm"]
  v$par3 <- 0
  v$par4 <- 0
  return(v)
}

#'chem.scenarios
#'
#' This function summarizes the all.scenarios data set by condensing each chemical-scenario combination into one row. The
#' number of rows of data for this combination on all.scenarios is recorded here.
#'
#' @param all This is the master list of all chemicals, and all exposure scenarios specific to each chemical, to be evaluated
#' in the current model run. This data set is compiled internally according to the user's specifications on the
#' \code{\link{Run_85687}} file. The user may create multiple scenarios files for special purposes, for example, for selected
#' chemical classes, or selected exposure pathways. A model run consists of two nested loops: the outer loop over chemicals and
#' the inner loop over the scenarios specific to that chemical.
#'
#' @return \item{chem.scen }{A data set consisting of all chemicals of interest, and the categories and scenarios for that chemical
#' with the potential for exposure.  The relationship between chemicals and scenarios is one to many. It is possible for one
#' chemical to have two or more "dermal" scenarios, provided each is in a different "category".  However, each combination of
#' chemical, category, and scenario must be unique. The last variable on chem.scen is "count", which is the number of rows of data
#' on all.scenarios devoted to this combination.}
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{Run}}, \code{\link{Source_chem_foods}}, \code{\link{Source_chem_prods}}, \code{\link{Source_scen_food}},  \code{\link{Source_scen_prods}}
#'
#' @keyword SHEDS
#'
#' @export

chem.scenarios = function(all) {
  chem.scen <- count(all,c("CAS","scenario","category"))
  names(chem.scen)[4] <- "count"
  chem.scen <- as.data.table(chem.scen)
  cat("\n Chemical-Scenario Table completed")
  return(chem.scen)
}

#' select.people
#'
#' Assigns demographic and physiological variables to each theoretical person to be modeled.
#'
#' @param n Number of persons.
#'
#' @param pop The population input; output of \code{\link{read.pop.file}} function. Contains counts by gender and each year of age
#' from the 2000 U.S. census.  When a large age range is modeled, this ensures that SHEDS chooses age and gender with the correct
#' overall probability.
#' #'
#' @param py Regression parameters on the three physiological variables of interest (weight, height, and body mass index) for
#' various age and gender groups. Output of the \code{\link{read.phys.file}} function.
#'
#' @param act.p Activity diary pools; output of the \code{\link{act.diary.pools}} function. Each element in this input is a list
#' of acceptable activity diary numbers for each year of age, for each gender, weekend, and season combination.
#'
#' @param diet.p Dietary diary pools; output of the \code{\link{diet.diary.pools}} function. Each element in this input is a list
#' of acceptable activity diary numbers for each year of age, for each gender, weekend, and season combination.
#'
#' @param act.d Activity diaries, which indicate the amount of time and level of metabolic activity in various 'micros'.
#' Each line of data represents one person-day (24 hours). Output of the \code{\link{read.act.diaries}} function.
#'
#' @param diet.d Daily diaries of dietary consumption by food group; output of the \code{\link{read.diet.diaries}} function.
#' Each line represents one person-day, with demographic variables followed by amounts (in grams/day) for a list of food types
#' indicated by a short abbreviation on the header line.
#'
#' @param specs Output of the \code{\link{read.run.file}} function, which can be modified by the \code{\link{update.specs}}
#' function before input into \code{\link{select.people}}.
#'
#' @details This is the first real step in the modeling process. It first fills an array \code{q} with uniform random numbers,
#' with ten columns because there are 10 random variables defined by this function. There number of rows correspond to the number
#' of persons, capped at the \code{set.size} specified in the \code{\link{Run}} input file (typically 5000).  Gender is selected
#' from a discrete (binomial) distribution where the counts of males and females in the study age range determines the gender
#' probabilities. Age is tabulated next, separately for each gender. The counts by year of age are chosen for the appropriate
#' gender and used as selection weights. Season is assigned randomly (equal weights) using those specified in the \code{\link{Run}}
#' input file. \code{Weekend} is set to one or zero, with a chance of 2/7 for the former.
#' The next block of code assigns physiological variables. Weight is lognormal in SHEDS, so a normal is sampled first and then
#' \code{exp()} is applied.  This means that the weight parameters refer to the properties of log(weight), which were fit by
#' linear regression. The basal metabolic rate (bmr), is calculated by regression. A minimum bmr is set to prevent extreme cases
#' from becoming zero or negative.  The alveolar breathing ventilation rate corresponding to bmr is also calculated. The SHEDS
#' logic sets activities in each micro to be a multiple of these rates, with outdoor rates higher than indoor, and indoor rates
#' higher than sleep rates. This calculated activities affect the inhaled dose. The skin surface area is calculated using
#' regressions based on height and weight for 3 age ranges.
#' The next step is to assign diaries.  Here, a FOR loop over n persons (n rows) is used to assign appropriate diet and activity
#' diary pools to each person. An empirical distribution is created,  consisting of the list of diary numbers for each pool.
#' The final step is to retrieve the actual data from the chosen activity and diet diaries, and the result becomes \code{pd}.
#'
#' @return pd A dataframe of "person-demographics": assigned demographic and physiological parameters for each theoretical
#' person modeled in SHEDS.HT
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{run}}, \code{\link{read.pop.file}}, \code{\link{read.phys.file}}, \code{\link{act.diary.pools}}, \code{\link{diet.diary.pools}}, \code{\link{read.act.diaries}}, \code{\link{read.diet.diaries}}, \code{\link{update.specs}}
#'
#' @keyword SHEDS
#'
#' @export

select.people = function(n, pop, py, act.p, diet.p, act.d, diet.d,specs) {
  # generate random numbers
  q <- matrix(runif(10*n),nrow=n,ncol=10)

  # find gender counts
  person  <- 1:n
  males   <- sum(pop$males)
  females <- sum(pop$females)
  if (males+females==0) stop("No population \n")
  gender <- distrib("disc",q=q[,1],v=c("F","M"),p=c(females,males))
  mode(pop$age) <- "character"
  age <- vector("numeric",n)
  if (females>0) {
    fage <- distrib("disc",q=q[,2],v=pop$age,p=pop$females)
    age[gender=="F"] <- fage[gender=="F"]
  }
  if (males>0)  {
    mage <- distrib("disc",q=q[,2],v=pop$age,p=pop$males)
    age[gender=="M"] <- mage[gender=="M"]
  }
  season  <- distrib("empi",q=q[,3],v=specs$seasons)
  weekend <- distrib("bern",q=q[,4],0.285714)

  # assign physiology
  r <- age+1+100*(gender=="M")
  py[is.na(py)]<- 0
  # weight is body weight in kg.
  weight <- exp(py$wgtmean[r]+py$wgtstdev[r]*distrib("norm",q=q[,5]))
  # height is in cm.
  ht1    <- py$hgtmean[r]+py$hgtstdev[r]*distrib("norm",q=q[,6])
  ht2    <- py$hgtinter[r]+py$hgtslope[r]*log(weight)+py$hgtresid[r]*
    distrib("norm",q=q[,6])
  height <- mapply(max,ht1,ht2)
  # bmr is basal metbolic rate in megajoules per day.
  bmr    <- py$bmrinter[r]+py$bmrslope[r]*weight+py$bmrresid[r]*
    distrib("norm",q=q[,7])
  bmr    <- mapply(max,0.1,bmr)
  # bva is basal ventilation rate in m3/day.
  bva    <- bmr*0.166*0.01963*1440*(0.20+0.01*q[,8])
  # surfarea is total skin surface area in cm2.
  surf1  <- 266.7 * height^0.38217 * weight^0.53937
  surf2  <- 305.0 * height^0.35129 * weight^0.54375
  surf3  <- 154.5 * height^0.54468 * weight^0.46366
  surfarea <- surf2
  surfarea[age<=5]  <- surf1[age<=5]
  surfarea[age>=20] <- surf3[age>=20]

  # select diary matched by gender, season, weekend, and age
  diary.id <- vector("integer",n)
  diet.id  <- vector("integer",n)
  for (i in 1:n) {
    atype <- paste0(gender[i],weekend[i],season[i],age[i])
    dtype <- paste0(gender[i],age[i])
    diary.id[i] <- distrib("empi",q=q[,9][[i]], v= act.p[[atype]])
    diet.id[i]  <- distrib("empi",q=q[,10][[i]],v=diet.p[[dtype]])
  }
  people <- as.data.table(data.frame(person,gender,age,season,weekend,
                                     weight,height,bmr,bva,surfarea,diary.id,diet.id))
  setkey(people,person,diary.id)
  pd1 <- join(people,act.d,by="diary.id")
  setkey(pd1,person,diet.id)
  pd  <- join(pd1,diet.d,by="diet.id")
  return(pd)
}

#' add.media
#'
#' For each theoretical person parameterized in the \code{pd} data frame, which is output from the \code{\link{select.people}}
#' function, this function generates the exposure duration for each potential exposure medium.
#'
#' @param n Number of persons.
#'
#' @param media List of the potential contact media for the model; output of the \code{\link{read.media.file}} function.
#' Each is found in a specific microenvironment (micro).
#'
#' @param pd Data frame containing internally assigned demographic and physiological parameters for each theoretical person
#' modeled in SHEDS.HT. Output of the \code{\link{select.people}} function
#'
#' @details To generate the exposure duration for each theoretical person, an array q of uniform random samples is generated with
#' rows = \code{n} (per set, where set size is specified by the user in the \code{\link{Run}} input file) and columns =
#' \code{nrow(media)} (the number of potential exposure media). An array called \code{dur} is created for the number of minutes on
#' the activity diary in the relevant micro multiplied by the relevant probability of contact (as specified in the
#' \code{media} input). The array consists of a row for each person and a column for each of the exposure media.  The output is a
#' data table containing the \code{pd} data frame and the \code{dur} array.
#'
#' @return pdm A data table containing the \code{pd} data frame of physiological and demographic parameters for each theoretical
#' person, and the \code{dur} array, which specifies the duration of exposure to each potential exposure medium for each person
#' in  \code{pd}.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{select.people}}, \code{\link{read.media.file}}
#'
#' @keyword SHEDS
#'
#' @export

add.media = function(n, media, pd) {
  # generate random numbers
  cols <- nrow(media)
  mic  <- data.frame(matrix(0,nrow=n,ncol=cols))
  dur  <- data.table(matrix(0,nrow=n,ncol=cols))
  setnames(dur,names(dur),paste0("dur.",tolower(media$media)))
  pdf  <- as.data.frame(pd)
  for (i in 1:cols){
    mic[,i]  <- pdf[,paste0(tolower(media$micro[i]),".min")]
    dur[,i]  <- mic[,i]*media$contact.p[i]
  }
  pdm <- data.table(pd,dur)
  return(pdm)
}

#' add.factors
#' Adds specific exposure factors to the \code{pdm} data table, which is output from the \code{\link{add.media}} function. Here,
#' "specific" means taking into account the age, gender, season, and exposure media for each person.
#'
#' @param n Number of persons
#'
#' @param gen.f Non-media specific exposure factors as a data table. Output from the \code{\link{gen.factor.tables}} function.
#'
#' @param med.f Media specific exposure factors presented as a data table. Output from the \code{\link{med.factor.tables}} function.
#'
#' @param exp.f Distributional parameters for the exposure factors. Output of the \code{\link{exp.factors}} function.
#'
#' @param surf A list of surface media. Modified output of the \code{\link{read.media.file}} function.
#'
#' @param pdm A data table containing the \code{pd} data frame of physiological and demographic parameters for each theoretical
#' person, and the \code{dur} array, which specifies the duration of exposure to each potential exposure medium for each person
#' in  \code{pd}. Output of the \code{\link{add.media}} function.
#'
#' @details The process of adding specific exposure factors to \code{pdm} involves multiple steps. First, \code{w} is determined,
#' which is the number of general factors plus the product of the number of media-specific factors and the number of surface media.
#' Air media do not have media specific factors in this version of SHEDS. An array, \code{q}, of uniform random samples is
#' generated, with one row per person and \code{w} columns.  A zero matrix, \code{r}, of the same size is defined.
#' Once these matrices are defined, the media-specific factors are determined. Two nested loops over variable and surface type
#' generate the values, which are stored in \code{r}. Next, another FOR loop determines the general factors. The \code{p} data set
#' contains the age, gender, and season for each person. These two data sets are then merged. The evaluation of these factors is
#' handled by the \code{\link{eval.factors}} function.
#' One of the exposure factors is \code{handwash.freq}. This was also part of SHEDS-Multimedia, where it represented the mean
#' number of hours in the day with hand washing events. An important aspect of that model was that because each person was followed
#' longitudinally, the actual number of hand washes on each day varied from one day to the next. Because of this, the distribution
#' for \code{handwash.freq} did not need to be restricted to integer values, as (for example) a mean of 4.5 per day is acceptable
#' and achievable, while choosing integer numbers of hand washes each day.  One of the early goals with SHEDS.HT was to attempt
#' to reproduce selected results from SHEDS-Multimedia. Therefore, similar logic was built into the current model. The
#' \code{hand.washes} variable is sampled from a distribution centered on \code{handwash.freq}, and then rounded to the nearest
#' integer.
#' The \code{bath} variable is another difficult concept. In theory, baths and showers are recorded on the activity diaries.
#' In practice, the activity diaries were constructed from approximately 20 separate studies, some of which did not contain enough
#' detail to identify separate bath or shower events. The result is that about half of all diaries record such events, but the
#' true rate in the population is higher. The \code{bath.p} variable was created to address this. It represents the  probability
#' that a non bath/shower activity diary should actually have one. Therefore, if the diary has one, then SHEDS automatically has
#' one. Otherwise, a binomial sample using \code{bath.p} as the probability is drawn. A bath/shower occurs unless both of these
#' are zero.
#' The effectiveness of hand washes or bath/shower at removing chemical from the skin is determined in the
#' \code{\link{post.exposure}} function.
#'
#' @return pdmf A data set containing the \code{pdm} data table as well as media specific exposure factors, the number of baths
#' taken, and the number of hand wash events occurring per day per person contained in \code{pdm}.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{eval.factors}}, \code{\link{post.exposure}}
#'
#' @keyword HEDS
#'
#' @export
#'
add.factors = function(n, gen.f, med.f, exp.f, surf, pdm) {
  cols <- length(surf)
  vmed <- as.character(unique(med.f$varname))
  vgen <- as.character(unique(gen.f$varname))
  w    <- length(vmed)*cols+length(vgen)
  med.f$media <- tolower(med.f$media)
  # generate random numbers
  q <- matrix(runif(n*w),nrow=n)
  r <- data.table(matrix(0,nrow=n,ncol=w))
  setkey(pdm,gender,season,age)
  for (v in 1:length(vmed)) {
    var <- vmed[v]
    for (m in 1:length(surf)) {
      med <- surf[m]
      vm <- med.f[med.f$varname==var & med.f$media==med]
      setkey(vm,gender,season)
      p <- join(pdm,vm,by=c("gender","season"))
      p <- p[(p$min.age<=p$age)&(p$age<=p$max.age)]
      setkey(p,person)
      col <- length(vmed)*(m-1)+v
      r[,col] <- eval.factors(p$row,q[,col],exp.f)
      setnames(r,paste0("V",col),paste0(var,".",med))
    }
  }
  for (v in 1:length(vgen)) {
    vg <- gen.f[gen.f$varname==vgen[v]]
    setkey(vg,gender,season)
    p <- join(pdm,vg,by=c("gender","season"))
    p <- p[(p$min.age<=p$age)&(p$age<=p$max.age)]
    setkey(p,person)
    col <- length(vmed)*length(surf)+v
    r[,col] <- eval.factors(p$row,q[,col],exp.f)
    setnames(r,paste("V",col,sep=""),vgen[v])
  }
  hand.washes <- round(r$handwash.freq +
                         sqrt(r$handwash.freq)*(2*q[,w]-1))
  bath <- as.logical(pdm$bath.mins+r$bath.p)
  pdmf <- cbind(pdm,r,bath,hand.washes)
  return(pdmf)
}

#' eval.factors
#'
#' Assigns distributions for each relevant exposure factor for each theoretical persons to be modeled in SHEDS.HT.
#'
#' @param r Vector of row numbers, equivalent to the number of theoretical people in the model sample. This input is created
#' internally by the \code{\link{add.factors}} function.
#'
#' @param q User-specified list of desired quantiles to be included in the output.
#'
#' @param ef Distributional parameters for the exposure factors; an output of the \code{\link{exp.factors}} function and an
#' input argument to the \code{\link{add.factors}} function
#'
#' @return z A data frame specifying the form of the distribution and relevant parameters for each combination of exposure
#' factor and individual person. The parameters include:
#' \describe{
#' \item{\code{form}}{The form of the distribution for each exposure factor (i.e., point, triangle).}
#' \item{\code{par1-par4}}{The parameters associated with the distributional form specified in form.}
#' \item{\code{lower.trun}}{the lower limit of values to be included in the distribution.}
#' \item{\code{upper.trun}}{The upper limit of values to be included in the distribution.}
#' \item{\code{resamp}}{Logical field where: yes=resample, no=stack at truncation bounds.}
#' \item{\code{q}}{The quantiles of the distribution corresponding to those specified by the user in the \code{q} input.}
#' \item{\code{p}}{The probabilities associated with the quantiles.}
#' \item{\code{v}}{The values associated with the quantiles.}
#' }
#'
#' @author Krisin Isaacs, Graham Glen
#'
#' @seealso \code{\link{add.factors}}, \code{\link{exp.factors}}
#'
#'@keyword SHEDS
#'
#' @export

eval.factors = function(r, q, ef) {
  n    <- length(r)
  ur   <- unique(r)
  nr   <- length(ur)
  z    <- vector("numeric",n)
  if (nr>0) { for (i in 1:nr) {
    j <- ur[i]
    z[r==j] <- distrib(ef$form[j],ef$par1[j],ef$par2[j],ef$par3[j],
                       ef$par4[j],ef$lower.trun[j],ef$upper.trun[j],
                       ef$resamp[j],q=q[r==j],p=ef$probs[j],v=ef$values[j])
  } }
  return(z)
}

#' make.cbase
#'
#' Extends the \code{base} data set, output from the \code{\link{generate.person.vars}} function, to include a set of
#' chemical-specific exposure variables. Currently, all new variables are initialized to zero.
#'
#' @param base Data set of all chemical-independent information needed for the exposure assessment.  Each row corresponds to
#' a simulated person. The variables consist of age, gender, weight, diet and activity diaries, food consumption, minutes in
#' each micro, and the evaluation of the exposure factors for that person.
#'
#' @param chem List of the chemical(s) of interest, determined via the \code{chemical} input in the \code{\link{Run}} file.
#' The exposure variables for the specified chemicals will be appended to the \code{base} data set.
#'
#' @return The output is a data set cb consisting of the base data set, appended by the following columns:
#'
#' \describe{
#' \item{chemrep }{The unique ID corresponding to each chemical specified in the \code{chem} argument.}
#' \item{inhal.abs.f }{For each person, the fraction of absorption of a given chemical when inhaled.}
#' \item{urine.f }{For each person, the fraction of intake of a given chemical that is excreted through urine.}
#' \item{exp.inhal.tot }{For each person, the total exposure of a given chemical via the direct inhalation scenario.}
#' \item{dose.inhal.tot }{For each person, the corresponding dose for total exposure of a given chemical via the direct
#' inhalation scenario.}
#' \item{exp.dermal.tot }{For each person, the  total exposure of a given chemical via the direct dermal application scenario.}
#' \item{exp.ingest.tot }{For each person, the  total exposure of a given chemical via the direct ingestion scenario.}
#' \item{exp.dietary.tot }{For each person, the total exposure of a given chemical via the dietary scenario.}
#' \item{exp.migrat.tot }{For each person, the total exposure of a given chemical via the migration scenario.}
#' \item{exp.nondiet.tot }{For each person, the total exposure of a given chemical via the non-dietary food exposure scenario.}
#' \item{exp.ddd.tot }{For each person, the total exposure of a given chemical via the down the drain scenario.}
#' \item{exp.window.tot }{For each person, the total exposure of a given chemical via the out the window scenario.}
#' }
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{Run}}, \code{\link{run}}
#'
#' @keyword SHEDS
#'
#' @export

make.cbase = function(base,chem) {
  n              <- nrow(base)
  chemrep        <- rep(chem,n)
  inhal.abs.f    <- rep(0.75,n)
  urine.f        <- rep(1,n)
  exp.inhal.tot  <- rep(0,n)
  exp.inhal.indir  <- rep(0,n)
  dose.inhal.tot <- rep(0,n)
  conc.inhal.max.prod.aer<- rep(0,n)
  conc.inhal.max.prod.vap<- rep(0,n)
  exp.dermal.tot <- rep(0,n)
  exp.ingest.tot <- rep(0,n)
  exp.dietary.tot<- rep(0,n)
  exp.migrat.tot <- rep(0,n)
  exp.nondiet.tot<- rep(0,n)
  exp.ddd.tot    <- rep(0,n)
  exp.window.tot <- rep(0,n)
  return(cbind(base,chemrep,inhal.abs.f,urine.f,exp.inhal.tot,
               dose.inhal.tot,exp.dermal.tot,exp.ingest.tot,exp.dietary.tot,
               exp.migrat.tot,exp.nondiet.tot,exp.ddd.tot,exp.window.tot,conc.inhal.max.prod.aer,conc.inhal.max.prod.vap,exp.inhal.indir))
}

#' create.scen.factors
#'
#' This function takes the information on distributions from the \code{all.scenarios} data set (which comes from
#' the \code{\link{source_variables_12112015}} input file) and converts it into the parameter set needed by SHEDS.HT.
#'
#' @param f An internally generated data set from the \code{\link{Source_vars}} input file (as specified in the \code{\link{Run}}
#' file) containing data on the distribution of each source variable in SHEDS.HT.
#'
#' @details The steps involved in this function are 1) converting \code{prevalence} from a percentage to a binomial form, 2)
#' converting \code{CV} to standard deviation for normals, and 3) converting \code{Mean} and \code{CV} to par1 and par2 for
#' lognormals.
#' Note that \code{prevalence} in SHEDS becomes a binomial distribution, which returns a value of either 0 or 1 when evaluated.
#' Each simulated person either "does" or "does not" partake in this scenario. Similar logic applies to the \code{frequency}
#' variable, except that the returned values may be larger than one (that is, 2 or more) for very frequent scenarios.  All the
#' exposure equations contain the \code{prevalence} variable. If \code{prevalence} is set to one for that person, the exposure
#' is as expected, but if \code{prevalence}=0 then no exposure occurs.
#'
#' @return dt A data set with values and probabilities associated with each variable distribution presented in \code{f}.
#' All variables in \code{f} are also retained.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{Source_vars}}, \code{\link{Run}}
#'
#' @keyword SHEDS
#'
#' @export

create.scen.factors = function(f) {
  options(warn=-1)
  mode(f$mean)    <- "numeric"
  mode(f$cv)      <- "numeric"
  mode(f$min.age) <- "numeric"
  mode(f$max.age) <- "numeric"
  options(warn=0)
  pars <- as.data.frame(matrix(NA,nrow=nrow(f),ncol=7))
  names(pars) <- c("par1","par2","par3","par4","lower.trun",
                   "upper.trun","row")
  pars$lower.trun <- 0
  pars$row        <- 1:nrow(pars)
  resamp          <- "yes"
  names(resamp)   <- "resamp"
  probs           <- ""
  names(probs)    <- "probs"
  if(!any(names(f)=="values")) f$values<-rep("",nrow(f))
  values          <- f$values
  names(values)   <- "values"
  for (i in 1:nrow(f)) {
    s <- strtrim(tolower(f$form[i]),4)
    if (tolower(f$varname[i])=="prevalence" &&
        tolower(f$form[i])=="point") {
      f$form[i]    <- "binomial"
      f$mean[i]    <- f$mean[i]/100
      pars$par1[i] <- f$mean[i]
    }
    if (s=="poin") { pars$par1[i] <- f$mean[i] }
    if (s=="norm") {
      pars$par1[i] <- f$mean[i]
      pars$par2[i] <- f$cv[i]*f$mean[i]
    }
    if (s=="logn") {
      cv <- as.numeric(f$cv[i])
      pars$par1[i] <- f$mean[i]/(1+cv^2)^0.5
      pars$par2[i] <- exp(log(1+cv^2)^0.5)
    }
  }
  vars <- c("varname","form","par1","par2","par3","par4",
            "lower.trun","upper.trun","resamp","min.age","max.age",
            "gender","row","values","probs")
  df <- as.data.frame(cbind(f,pars,resamp,values,probs))[vars]
  dt <- as.data.table(df)
  return(dt)
}

#' scen.factor.indices
#'
#' Constructs scenario-specific exposure factors for each relevant combination of age and  gender. The input is derived
#' internally from the \code{\link{Source_vars}} file specified on the \code{\link{Run}} file.
#'
#' @param sdat The chemical-scenario data specific to a given combination of chemical and scenario.
#'
#' @param expgen Non-media specific exposure factors as a data table. Output from \code{\link{gen.factor.tables}}
#'
#' @details The constructed factors are not media-specific, although they are scenario-specific,  and most scenarios include
#' just one surface medium.
#' This function evaluats the scenario-specific exposure factors separately from the other factors because these
#' scenario-specific factors may change with every chemical and scenario, whereas the other factors remain the same across
#' chemicals and scenarios for each person.
#'
#' @return Returns indices of scenario-specific exposure factors for each relevant age and gender combination. The output
#' is only generated internally.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @keyword SHEDS
#'
#' @export

scen.factor.indices = function(sdat,expgen) {
  vars <- c("varname","min.age","max.age","gender")
  dt <- as.data.table(as.data.frame(sdat)[vars])
  rgen <- nrow(dt)
  dt$row <- 1:rgen
  expgen <- dt
  for (i in 1:nrow(dt)) {
    if (is.na(dt$gender[i])|dt$gender[i]=="") { g <- c("F","M")
    } else g <- as.character(dt$gender[i])
    for (ig in 1:length(g)) {
      rgen <- rgen+1
      expgen <- rbind(expgen,dt[i])
      expgen$gender[[rgen]] <- g[[ig]]
    }
  }
  return(expgen[(nrow(dt)+1):nrow(expgen)])
}

#' dir.dermal
#'
#' Models the dermal exposure scenario for each theoretical person.
#'
#' @param sd The chemical-scenario data specific to relevant combinations of chemical and scenario. Generated internally.
#'
#' @param cd The list of scenario-specific information for the chemicals being evaluated. Generated internally.
#'
#' @details The dermal exposure scenario is relatively straightforward. The function produces a prevalence value, which reflects
#' the fraction of the population who use this scenario at all. It also produces a frequency value, which is the mean number of
#' times per year this scenario occurs among that fraction of the population specified by prevalence.
#' Since SHEDS operates on the basis of one random day, the frequency is  divided by 365 and then passed to the
#' \code{\link{p.round}} (probabilistic rounding) function, which rounds either up or down to the nearest integer. Very common
#' events may happen more than once in a day.
#' The function also produces a mass variable, which refers to the  mass of the product in grams in a typical usage event.
#' The composition is the percentage of that mass that is the chemical in question.
#' The resid variable measures the fraction that is likely to remain on the skin when the usage event ends.
#' The final output is the total dermal exposure for each chemical-individual combination in micrograms, which is the product of
#' the above variables multiplied by a factor of 1E+06.
#'
#' @return dir.derm  For each person, the calculated direct dermal exposure that occurs when a chemical-containing product
#' is used. This does not include later contact with treated objects, which is indirect exposure.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{run}}, \code{\link{p.round}}
#'
#' @keyword SHEDS
#'
#' @export

dir.dermal = function(sd,cd) {
  n <- nrow(sd)
  exp.inhal.dirderm  <- rep(0,n)
  dose.inhal.dirderm <- rep(0,n)
  exp.ingest.dirderm <- rep(0,n)
  if (exists("f.contact",sd)==FALSE) sd$f.contact <- 1
  exp.dermal.dirderm <- sd$use.prev * sd$use.today *
    sd$mass * cd$chem.prev * cd$f.chemical *
    sd$f.residual * sd$f.contact * 1E6
  # mass in [g], use.freq/365 gives #/day, others are [-], 1E6 is [ug/g]
  dir.derm  <- as.data.table(cbind(exp.dermal.dirderm, exp.ingest.dirderm,
                                   exp.inhal.dirderm, dose.inhal.dirderm))
  return(dir.derm)
}

#' dir.ingested
#'
#' Models the ingestion exposure scenario for each theoretical person.
#'
#' @param sd The chemical-scenario data specific to relevant combinations of chemical and scenario. Generated internally.
#'
#' @param cd The list of scenario-specific information for the chemicals being evaluated. Generated internally.
#'
#' @details This scenario is for accidental ingestion during product usage. Typical examples are toothpaste, mouthwash,
#' lipstick or chap stick, and similar products used on the face or mouth.
#' The function produces a \code{prevalence} value, which reflects the fraction of the population who use this scenario at all.
#' It also produces a \code{frequency} value, which is the mean number  of times per year this scenario occurs among that
#' fraction of the population specified by prevalence.
#' Since SHEDS operates on the basis of one random day, the frequency is  divided by 365 and then passed to the
#' \code{\link{p.round}} (probabilistic rounding) function, which rounds either up or down to the nearest integer. Very common
#' events may happen more than once in a day.
#' The function also produces a \code{mass} variable, which refers to the  mass of the product in grams in a typical usage event.
#' The \code{composition} is the percentage of that mass that is the chemical in question.
#' The \code{ingested} variable represents the percentage of the mass applied that becomes ingested. Since these products are
#' not intended to be swallowed, this should typically be quite small (under 5\%).
#' The final output is the total incidental ingested exposure for each chemical-individual combination in micrograms, which
#' is the product of the above variables multiplied by a factor of 1E6.
#'
#' @return dir.ingest For each person, the calculated quantity of a given chemical incidentally ingested during or
#' immediately after use of products such as toothpaste. Does not include exposure via food and drinking water.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @keyword SHEDS
#'
#' @export

dir.ingested = function(sd,cd) {
  n <- nrow(sd)
  exp.inhal.diringest  <- rep(0,n)
  dose.inhal.diringest <- rep(0,n)
  exp.dermal.diringest <- rep(0,n)
  exp.ingest.diringest <- sd$use.prev * sd$use.today *
    sd$mass * cd$chem.prev * cd$f.chemical *
    sd$f.ingested * 1E6
  # mass in [g], chem.prev & f.chemical & f.ingested are [-], 1E6 gives [ug]
  dir.ingest  <- as.data.table(cbind(exp.dermal.diringest,
                                     exp.ingest.diringest,exp.inhal.diringest, dose.inhal.diringest))
  # cat("\n dir.ingest = ",dir.ingest$exp.ingest.diringest[1:n])
  return(dir.ingest)
}

#' dir.inhal.aer
#'
#' Models the inhalation exposure from the use of aerosol products for each theoretical person.
#'
#' @param sd The chemical-scenario data specific to relevant combinations of chemical and scenario. Generated internally.
#'
#' @param cd The list of scenario-specific information for the chemicals being evaluated. Generated internally.
#'
#' @param cb A copy of the \code{base} data set output from the \code{\link{make.cbase}} function, with columns added for
#' exposure variables.
#'
#' @param io A binary variable indicating whether the volume of the aerosol is used to approximate the affected volume.
#'
#' @details This scenario considers inhalation exposure from the use of aerosol products. Typical examples include hairspray
#' and spray-on mosquito repellent.
#' The function produces a \code{prevalence} value, which reflects the fraction of the population who use this scenario at all.
#' It also produces a \code{frequency} value, which is the mean number  of times per year this scenario occurs among that
#' fraction of the population specified by prevalence.
#' Since SHEDS operates on the basis of one random day, the \code{frequency} is  divided by 365 and then passed to
#' the \code{\link{p.round}} (probabilistic rounding) function, which rounds either up or down to the nearest integer. Very
#' common events may happen more than once in a day.
#' The function also produces a \code{mass} variable, which refers to the  mass of the product in grams in a typical usage
#' event. The \code{composition} is the percentage of that mass that is the chemical in question.
#' \code{frac.aer} is the fraction of the product mass that becomes aerosolized, and the \code{volume} affected by the use is
#' approximated to allow the calculation of a concentration or density. Defaults are set in the code if these variables  are
#' missing from the input file.
#' Exposure for the inhalation pathway has units of micrograms per cubic meter, reflecting the average air concentration of
#' the chemical. An \code{airconc} variable is defined using \code{mass}, \code{composition}, \code{frac.aer}, and \code{volume}.
#' Since exposure depends on the time-averaged concentration, a duration is necessary. For example, if one spends five minutes
#' in an aerosol cloud and the rest of the day in clean air, the daily exposure is the cloud concentration multiplied by 5/1440
#' (where 1440 is the number of minutes in a day).
#' This function also calculates the inhaled dose, in units of micrograms per day. The dose equals the product of
#' \code{exposure} (g/m3), basal ventilation rate, \code{bvr} (m3/day), the METS factor of 1.75 (typically people inhale air
#' at an average of 1.75 times the basal rate to support common daily activities), and a  conversion factor of 1E+06 from grams
#' to micrograms.
#'
#' @return dir.inh.aer The calculated quantity of chemical inhalation from aerosols, like hairspray and similar products,
#' that are directly injected into the air on or around each exposed person.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' keyword{ ~SHEDS }
#'
#' @export

dir.inhal.aer = function(sd,cd,cb,io) {
  n <- nrow(sd)
  exp.dermal.dirinhaer <- rep(0,n)
  exp.ingest.dirinhaer <- rep(0,n)
  if (exists("f.aerosol",sd)==FALSE) sd$f.aerosol <- 0.5
  if (exists("volume"  ,sd)==FALSE && io==1)  sd$volume <- 24
  if (exists("volume"  ,sd)==FALSE && io==0)  sd$volume <- 480
  if (exists("duration",sd)==FALSE)  sd$duration  <- 1
  volume               <- sd$volume
  if (io==1) volume[volume<24]  <- 24
  if (io==0) volume[volume<480] <- 480
  airconc <- sd$mass * cd$chem.prev * cd$f.chemical * sd$f.aerosol / volume
  # the following limits chemical concentrations to 0.1% of air
  airconc[airconc>1.2] <- 1.2
  # mass in g, vol in [m3], others unitless, airconc in [g/m3]
  exp.inhal.dirinhaer  <- sd$use.prev * sd$use.today *
    airconc * sd$duration/1440 *1E6
  conc<-sd$use.prev * sd$use.today *airconc* 1E6
  # exp.inhal in [ug/m3], dose.inhal in [ug/day]
  #dose.inhal.dirinhaer <- exp.inhal.dirinhaer*cb$bva*1.75
  dose.inhal.dirinhaer <- exp.inhal.dirinhaer*cb$bva*cb$in.awk.pai
  # Typical breathing rates are about 1.75 * basal rate in [m3/day]
  # Changed by WGG on Apr 13, 2016 to PAI
  dir.inh.aer <- as.data.table(cbind(exp.dermal.dirinhaer,
                                     exp.ingest.dirinhaer, exp.inhal.dirinhaer, dose.inhal.dirinhaer,conc))
  if(min(exp.inhal.dirinhaer)<0) cat("\n Negative exposure dir.inhal.aer")
  return(dir.inh.aer)
}

#' dir.inhal.vap
#'
#' Models the inhalation exposure from the vapors of volatile chemicals for each theoretical person.
#'
#' @param sd The chemical-scenario data specific to relevant combinations of chemical and scenario. Generated internally.
#'
#' @param cd The list of scenario-specific information for the chemicals being evaluated. Generated internally.
#'
#' @param cprops The chemical properties required for SHEDS-HT. The default file (the \code{\link{Chem_props file}} read in by
#' the \code{\link{read.chem.props}} function and modified before input into the current function) was prepared from publicly
#' available databases using a custom program (not part of SHEDS-HT).  The default file contains 7 numerical inputs per
#' chemical, and the required properties are molecular weight (\code{MW}), vapor pressure (\code{VP.Pa}), solubility
#' (\code{water.sol.mg.l}), octanol-water partition coefficient (\code{log.Kow}), air decay rate (\code{half.air.hr}),
#' decay rate on surfaces (\code{half.sediment.hr}), and permeability coefficient (\code{Kp}).
#'
#' @param cb A copy of the \code{base} data set output from the \code{\link{make.cbase}} function, with columns added
#' for exposure variables.
#'
#' @param io A binary variable indicating whether the volume of the aerosol is used to approximate the affected volume.
#'
#' @details This scenario considers inhalation exposure from vapors (not aerosols). For example, painting will result in the
#' inhalation of vapor, but it does not involve aerosols (unless it is spray paint).
#' For this scenario, the vapor pressure and the molecular weight are relevant variables for determining exposure. These
#' variables are included in the input to the \code{cprops} argument, which is drawn internally from the
#' \code{\link{Chem_props}} file.
#' The function produces a \code{prevalence} value, which reflects the fraction of the population who use this scenario at
#' all. It also produces a \code{frequency} value, which is the mean number of times per year this scenario occurs among that
#' fraction of the population specified by \code{prevalence}.
#' Since SHEDS operates on the basis of one random day, the \code{frequency} is  divided by 365 and then passed to the
#' \code{\link{p.round}} (probabilistic rounding) function, which rounds either up or down to the nearest integer. Very
#' common events may happen more than once in a day.
#' The function also produces a \code{mass} variable, which refers to the  mass of the product in grams in a typical usage
#' event. The \code{composition} is the percentage of that mass that is the chemical in question. The \code{evap} variable is
#' an effective evaporated mass, calculated using the \code{mass}, \code{composition} (converted from percent to a fraction),
#' the vapor pressure as a surrogate for partial pressure, and \code{duration} of product use.  The \code{duration} term is
#' made unitless by dividing by 5 (minutes), which is an assumed time constant.
#' The effective air concentation \code{airconc} is calculated as \code{evap}/\code{volume}. The value for \code{airconce}
#' is capped by \code{maxconc}, which represents the point at which evaporation ceases.  For chemicals used for a short
#' duration, or with low vapor presssure, \code{maxconc} might not be reached before usage stops.
#' Once \code{airconc} is established, the function also calculates the inhaled dose, in units of micrograms per day.
#' The dose equals the product of \code{exposure} (g/m3), basal ventilation rate, \code{bvr} (m3/day), the METS factor of 1.75
#' (typically people inhale air at an average of 1.75 times the basal rate to support common daily activities), and a
#' conversion factor of 1E6 from grams to micrograms.
#'
#' @retrun dir.inh.vap The calculated quantity of chemical inhalation from exposure to vapors, such as those emitted by paint.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{run}}, \code{\link{p.round}}, \code{\link{read.chem.props}}, \code{\link{Chem_props}}
#'
#' @keyword SHEDS
#'
#' @export

dir.inhal.vap = function(sd,cd,cprops,cb,io) {
  n <- nrow(sd)
  exp.dermal.dirinhvap <- rep(0,n)
  exp.ingest.dirinhvap <- rep(0,n)
  if (exists("volume"  ,sd)==FALSE && io==1)  sd$volume <- 24
  if (exists("volume"  ,sd)==FALSE && io==0)  sd$volume <- 480
  vapor                <- cprops$vapor/1E5
  # vapor is converted to partial pressure
  molwt                <- cprops$mw
  if (is.na(vapor))  vapor <- 1
  if (is.na(molwt))  molwt <- 100
  vapfrac              <- vapor
  # vapfrac is the mass fraction that can evaporate
  vapfrac[vapfrac>1]   <- 1
  evap                 <- sd$mass * cd$chem.prev * cd$f.chemical * vapfrac
  # mass in [g], evap in [g]
  volume               <- sd$volume
  if (io==1) volume[volume<24]  <- 24
  if (io==0) volume[volume<480] <- 480
  airconc              <- evap/volume
  # the following limits chemical concentrations to 0.1% of air
  airconc[airconc>1.2] <- 1.2
  exp.inhal.dirinhvap  <- sd$use.prev * sd$use.today *
    airconc * (sd$duration/1440) * 1E6
  conc<-sd$use.prev * sd$use.today *airconc* 1E6
  # exp.inhal in [ug/m3], dose.inhal in [ug/day]
  #dose.inhal.dirinhvap <- exp.inhal.dirinhvap*cb$bva*1.75
  dose.inhal.dirinhvap <- exp.inhal.dirinhvap*cb$bva*cb$in.awk.pai
  # bva in [m3/day], 1.75 is PAI [-], dose in [ug/day]
  # Changed by WGG on Apr 13, 2016 to PAI
  dir.inh.vap <- as.data.table(cbind(exp.dermal.dirinhvap,
                                     exp.ingest.dirinhvap,exp.inhal.dirinhvap, dose.inhal.dirinhvap,conc))
  if(min(exp.inhal.dirinhvap)<0) cat("\n Negative exposure Inhalvap")
  return(dir.inh.vap)
}

#' down.the.drain.mass
#'
#' Models the quantity of chemical entering the waste water system on a per person-day basis.
#' @param sd The chemical-scenario data specific to relevant combinations of chemical and scenario. Generated internally.
#'
#' @param cd The list of scenario-specific information for the chemicals being evaluated. Generated internally.
#'
#' @details This function models the simplest of all the current scenarios. It evaluates the amount of chemical entering
#' the waste water system, on a per person-day basis.  The "exposure" is to a system, not a person, but this method uses
#' one person's actions to estimate their contribution to the total.
#' The function produces a \code{prevalence} value, which reflects the fraction of the population who use this scenario at
#' all. It also produces a \code{frequency} value, which is the mean number  of times per year this scenario occurs among
#' that fraction of the population specified by prevalence.
#' Since SHEDS operates on the basis of one random day, the \code{frequency} is  divided by 365 and then passed to the
#' \code{\link{p.round}} (probabilistic rounding) function, which rounds either up or down to the nearest integer. Very
#' common events may happen more than once in a day.
#' The function also produces a \code{mass} variable, which refers to the  mass of the product in grams in a typical
#' usage event. The \code{composition} is the percentage of that mass that is the chemical in question.
#' The final output, \code{exp.ddd.mass}, is the product of the \code{prevalence}, \code{frequency}, \code{mass},
#' \code{composition}, and the fraction going down the drain (\code{f.drain}, a variable in the \code{sd} input).
#' The result is in grams per person-day.
#'
#' @return exp.ddd.mass The calculated quantity of chemical going down the drain (i.e., from laundry detergent) and
#' entering the sewer system per person per day.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{run}}, \code{\link{p.round}}
#'
#' @keyword SHEDS
#'
#' @export

down.the.drain.mass = function(sd,cd) {
  n <- nrow(sd)
  exp.ddd.mass <- rep(0,n)
  exp.ddd.mass <- sd$use.prev * sd$use.today *
    sd$mass * cd$chem.prev * cd$f.chemical * sd$f.drain
  # mass in g, others unitless, ddd.mass in [g/day]
  return(signif(exp.ddd.mass,5))
}

#' food.residue
#'
#' Models exposure to chemicals from consumption of food containing a known chemical residue for each theoretical person.
#'
#' @param cdata The list of scenario-specific information for the chemicals being evaluated. Generated internally.
#'
#' @param cb Output of the \code{\link{make.cbase}} function.
#'
#' @param ftype Food consumption database which stores data on consumption in grams per day of each food type for each
#' person being modeled. Generated internally from the  \code{\link{Diet_diaries}} data set.
#'
#' @details In this function, the variable \code{foods} is defined as a list of the names of the food groups,
#' which are stored in both the \code{ftype} and \code{cb} arguments. The FOR loop picks the name of each food group as
#' a string, and converts it to a variable name. For example, the food group "FV" may have corresponding variables
#' \code{residue.FV}, \code{zeros.FV} (number of nondetects), and \code{nonzeros.FV} (number of detects) on the \code{cb}
#' data set. If these variables are not present, no  exposure results.  The variables may be present, but an individual
#' person may still receive zero exposure because  not all samples of that food are contaminated, as indicated by the
#' \code{zeros.FV} value.  The \code{nonzeros.FV} value is used to  determine the likelihood of contamination, with the
#' residue variable determining the amount found (when it is nonzero).
#' The dietary exposure is the product of the \code{consumption} as reported in the \code{ftypes} input (in grams) and
#' the \code{residue.FV} (in micrograms of chemical per gram of food), summed over all food groups. Three other vectors are
#' returned (corresponding to dermal exposure, inhalation exposure, and inhalation dose), but these are currently set to zero
#' for the dietary pathway. In principle, eating food could result in dermal exposure (for finger foods), or inhalation (as
#' foods may have noticeable odors), but such exposures are not large enough to be of concern at present.
#'
#' @return dietary The calculated quantity of chemical exposure from food residue in grams per person per day.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{Diet_diaries}}, \code{\link{run}}, \code{\link{p.round}}, \code{\link{generate.person.vars}}, \code{\link{food.migration}}
#'
#' @keyword SHEDS  kwd2
#'
#' @export

food.residue = function(cdata,cb,ftype) {
  n <- nrow(cb)
  exp.inhal.dietary  <- rep(0,n)
  dose.inhal.dietary <- rep(0,n)
  exp.dermal.dietary <- rep(0,n)
  diet               <- rep(0,n)
  names(diet)        <- "exp.diet"
  if (exists("residue",cdata)==TRUE) {
    q    <- runif(n)
    cbf  <- eval(parse(text=paste0("cb$",tolower(ftype))))
    if (!is.null(cbf)) {
      frc  <- cdata$detects/(cdata$detects+cdata$nondetects)
      diet <- diet + ifelse(q<frc, cbf*cdata$residue, 0)
    }
  }
  # residue in [ug/g], food mass cbf in [g], result in [ug]
  exp.ingest.dietary <- diet
  if(min(exp.ingest.dietary)<0) cat("\n Negative exposure dir.dietary")
  dietary  <- as.data.table(cbind(exp.dermal.dietary, exp.ingest.dietary,
                                  exp.inhal.dietary, dose.inhal.dietary))
  return(dietary)
}

#' food.migration
#'
#' Calculates chemical exposure from migration of chemicals from packaging or other contact materials into food (which is
#' then consumed). This is added to the \code{dietary} output calculated by the \code{\link{food.residue}} function to
#' establish a total exposure from the dietary pathway.
#'
#' @param cdata The list of scenario-specific information for the chemicals being evaluated. Generated internally.
#'
#' @param sdata The chemical-scenario data specific to relevant combinations of chemical and scenario. Generated internally.
#'
#' @param cb A copy of the \code{base} data set output from the \code{\link{make.cbase}} function, with columns added for
#' exposure variables.
#'
#' @param ftype Food consumption database which stores data on consumption in grams per day of each food type for each person
#' being modeled. Generated internally from the  \code{\link{Diet_diaries}} data set.
#'
#' @details If migration data is stored in the \code{cdata} argument, this is added to the total exposure calculated in the
#' \code{\link{food.residue}} function.
#'
#' @return dietary The calculated quantity of chemical exposure from food residue and migration of chemicals into food from
#' packaging and other contact materials in grams per person per day.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{Diet_diaries}}, \code{\link{run}}, \code{\link{p.round}}, \code{\link{make.cbase}}, \code{\link{food.residue}}
#'
#' @keyword SHEDS
#'
#' @export

food.migration = function(cdata,sdata,cb,ftype) {
  n <- nrow(cb)
  exp.inhal.migrat  <- rep(0,n)
  dose.inhal.migrat <- rep(0,n)
  exp.dermal.migrat <- rep(0,n)
  migrat            <- rep(0,n)
  chem.prev         <- rep(1,n)
  contact           <- rep(1,n)
  # if chem.prev and contact not on input files, set to 1 by default
  if(exists("chem.prev",cdata)) chem.prev <- cdata$chem.prev
  if(exists("contact",sdata))   contact   <- sdata$contact
  if (exists("migration.conc",cdata)==TRUE) {
    q    <- runif(n)
    cbf  <- eval(parse(text=paste0("cb$",tolower(ftype))))
    if (!is.null(cbf)) {
      migrat <- migrat + cbf*cdata$migration.conc*p.round(chem.prev)*contact
    }
  }
  # migration [ug/g], food mass cbf [g], chem.prev [-], contact [-], result [ug]
  exp.ingest.migrat <- migrat
  if(min(exp.ingest.migrat)<0) cat("\n Negative exposure food.migration")
  dietary  <- as.data.table(cbind(exp.dermal.migrat, exp.ingest.migrat,
                                  exp.inhal.migrat, dose.inhal.migrat))
  return(dietary)
}

#' indir.exposure
#'
#'Models the indirect exposure to chemicals in the home for each theoretical person.
#'
#' @param sd The chemical-scenario data specific to relevant combinations of chemical and scenario. Generated internally.
#'
#' @param cb A copy of the \code{base} data set output from the \code{\link{make.cbase}} function, with columns added for
#' exposure variables.
#'
#' @param concs The concentration of the chemical (in air and/or on surfaces) being released into the environment. Outpu of
#' the \code{\link{get.fug.concs}} function.
#'
#' @param chem.data The list of scenario-specific information for the chemicals being evaluated. Generated internally.
#'
#' @details Indirect exposure happens after a product is no longer being used or applied, due to chemical lingering on various
#' surfaces or in the air. People who come along later may receive dermal or inhalation exposure from residual chemical in the
#' environment.
#' SHEDS.HT currently has two indirect exposure scenarios. One which applies to a one-time chemical treatment applied to a
#' house, and another which applies to continual releases from articles. Both scenarios consist of two parts: the first
#' determines the appropriate air and surface concentrations. That code is in the \code{\link{Fugacity}} module.  The second
#' part is the exposure calculation by the current function. Both types of indirect exposure scenarios call this function.
#' The surface and air concentrations from the \code{\link{Fugacity}} module are premised on the product use actually occurring.
#' Hence, \code{\link{indir.exposure}} starts by multiplying those concentrations by the \code{prevalence} (which is either 0
#' or 1, evaluated separately for each person).
#' For air, the exposure is the average daily concentration, which is the event concentration multiplied by the fraction of
#' the day spent in that event. The inhaled dose is the product of the exposure, the basal ventilation rate, and the PAI factor
#' (multiplier for the basal rate). A factor of 1E+06 converts the result from grams per day to micrograms per day.
#' Dermal exposure results from skin contact with surfaces. The surface concentration (ug/cm2) is multiplied by the
#' fraction available for transfer (\code{avail.f}, unitless), the transfer coefficient (\code{dermal.tc}, in cm2/hr), and the
#' contact \code{duration} (hr/day).  The result is the amount of chemical transferred onto the skin (ug/day).
#'
#' @return indir Indirect exposure to chemicals in the home in ug per person per day.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso \code{\link{Fugacity}}, \code{\link{get.fug.concs}}, \code{\link{make.cbase}}, \code{\link{run}}
#'
#' @keyword SHEDS
#'
#' @export

indir.exposure = function(sd,cb,concs,chem.data) {
  if (exists("home.prev",sd)) home.prev<-sd$home.prev else home.prev <- sd$use.prev
  surface             <- home.prev * chem.data$chem.prev * concs$conc.sur    # surface [ug/m2]
  air                 <- home.prev * chem.data$chem.prev * concs$conc.air    # air [ug/m3]
  exp.inhal.indirect  <- air * (cb$dur.in.awk.air+cb$dur.in.slp.air)/1440    # exp.inhal,indirect [ug/m3]
  dose.inhal.indirect <- exp.inhal.indirect * cb$bva * 1.75                  # dose.inhal.indirect [ug/day]
  dosefactor          <- (cb$dur.in.awk.air*cb$in.awk.pai+cb$dur.in.slp.air*cb$in.slp.pai)/1440
  dose.inhal.indirect <- cb$bva * air * dosefactor                           # dose.inhal.indirect [ug/day]
  exp.dermal.indirect <- surface * cb$avail.f.in.awk.srf /10000 *
    cb$dermal.tc.in.awk.srf/60 * cb$dur.in.awk.srf      # exp.dermal.indirect [ug/day]
  # /60 converts tc [cm2/hr] to [cm2/min], /10000 converts [cm2] to [m2], dur in [min/day]
  exp.ingest.indirect <- exp.dermal.indirect * cb$om.ratio.in.awk.srf        # exp.ingest.indirect [ug/day]
  indir <- as.data.table(cbind(exp.dermal.indirect,exp.ingest.indirect,
                               exp.inhal.indirect,dose.inhal.indirect))
  return(indir)
}

#' post.exposure
#'
#'Resolves the fate of chemical after initial exposure has occurred. This involves parsing out the amount of chemical
#'removed and amount ultimately contributing to the exposure dose for each person.
#'
#' @param cb A copy of the \code{base} data set output from the \code{\link{make.cbase}} function, with columns added for
#' exposure variables.
#'
#' @param cprops The chemical properties required for SHEDS-HT. The default file (the \code{\link{Chem_props file}} read in
#' by the \code{\link{read.chem.props}} function and modified before input into the current function) was prepared from
#' publicly available databases using a custom program (not part of SHEDS-HT).  The default file contains 7 numerical inputs
#' per chemical, and the required properties are molecular weight (\code{MW}), vapor pressure (\code{VP.Pa}), solubility
#' (\code{water.sol.mg.l}), octanol-water partition coefficient (\code{log.Kow}), air decay rate (\code{half.air.hr}), decay
#' rate on surfaces (\code{half.sediment.hr}), and permeability coefficient (\code{Kp}).
#'
#' @details This function resolves the fate of chemical after initial exposure has occurred.  The same function applies to all
#' exposure scenarios.
#' The dermal exposure is the most complicated, as there are five removal methods. All five are randomly sampled. While the sum
#' of the five means is close to one, the sum of five random samples might not be, so these samples are treated as fractions of
#' their sum. The \code{rem.bath} variable is either 0 or 1, so the bath removal term is either zero or the sampled bath removal
#' efficiency. The \code{rem.brush} term is also simple. The handwashing (\code{rem.wash}) and hand-to-mouth
#' (\code{rem.hmouth}) transfer terms are non-linear, because higher frequencies have less chemical available for removal on
#' each repetition.  The algorithms in place were fitted to output from SHEDS-Multimedia, which were summed to daily totals.
#' The final removal term is dermal absorption (\code{rem.absorb}). The base value is multiplied by the \code{Kp} factor from
#' the \code{cprops} argument, and divided by the value for permethrin, as the values from SHEDS-Multimedia were based on a
#' permethrin run.
#' The five terms are evaluated separately for each person, as is their sum. Each is then converted to a fraction of the whole.
#' The fractions may be quite different from one person to another. For one, perhaps 80\% of the dermal loading is removed by
#' a bath/shower, while for another person it is 0\% because they did not take one.  For the latter person, other four removal
#' terms are (on average) five times larger than for the former person, because together they account for 100\% of the removal,
#' instead of just 20\%.
#' The rest of the \code{post.exposure} function is mostly a matter of bookkeeping. The hand-to-mouth dermal removal term
#' becomes an ingestion exposure term. Summing exposures across routes is dubious, in part because inhalation exposures use
#' different units from the others, and because much of the dermal exposure never enters the body. In addition, summing dermal
#' and ingestion exposures may double-count the hand-to-mouth term. However, intake dose may be summed. In SHEDS, "intake dose"
#' is the sum of the inhaled dose (which is the amount of chemical entering the lungs in ug/day), the ingestion exposure (which
#' is the amount entering the GI tract in ug/day), and the dermal absorption (the amount penetrating into or through the skin,
#' so it cannot otherwise be removed, in ug/day).
#' The "absorbed dose" is also calculated: for dermal it is the same as the intake dose, but for ingestion and
#' inhalation there is another absorption factor, which was set on the exposure factors input file. An estimate of the
#' chemical in urine (in ug/day) is made.  Both the intake dose and the absorbed dose are reported in both (ug/day) and in
#' (mg/kg/day). Note that the latter requires the body weights of each individual. These cannot be obtained from the former
#' just by knowing the average body weight in SHEDS.
#' The above variables for each simulated person are written to the \code{fexp} object (the name stands for "final exposure").
#' Each chemical writes over the previous \code{fexp}, so the data must first be summarized and written to an output file.
#'
#' @return fexp Absorbed dose and intake dose of a given chemical for each theoretical person being modeled for all
#' exposure scenarios.
#'
#' @author Kristin Isaacs, Graham Glen
#'
#' @seealso  \code{\link{make.cbase}}, \code{\link{Chemprops_small}}
#'
#' @keyword SHEDS kwd2
#'
#' @export
#'
post.exposure = function(cb,cprops) {
  remf.bath   <- cb$rem.bath.f * cb$bath
  remf.wash   <- cb$rem.handwash.f * (1-exp(-cb$hand.washes/1.6))
  remf.hmouth <- cb$rem.handwash.f * exp(-cb$hand.washes/11) *
    cb$handmouth.area.f * sqrt(cb$handmouth.freq)
  remf.brush  <- cb$rem.brushoff.f
  Kpfactor    <- cprops$kp/0.208
  if(is.na(Kpfactor)) Kpfactor <- 1
  # default Kp is value for Permethrin
  remf.absorb <- cb$rem.derm.abs.f * Kpfactor
  sum <- remf.bath + remf.wash + remf.hmouth + remf.brush + remf.absorb

  rem.bath   <- cb$exp.dermal.tot * remf.bath/sum
  rem.wash   <- cb$exp.dermal.tot * remf.wash/sum
  rem.hmouth <- cb$exp.dermal.tot * remf.hmouth/sum
  rem.brush  <- cb$exp.dermal.tot * remf.brush/sum
  rem.absorb <- cb$exp.dermal.tot * remf.absorb/sum

  ingest.abs.f<-cprops$fabs
  exp.hmouth.tot   <- rem.hmouth
  exp.ingest.tot   <- cb$exp.ingest.tot + exp.hmouth.tot
  dose.intake.ug   <- exp.ingest.tot + cb$dose.inhal.tot + rem.absorb
  dose.intake.mgkg <- dose.intake.ug/(1000*cb$weight)
  abs.dermal.ug    <- rem.absorb
  abs.ingest.ug    <- exp.ingest.tot * ingest.abs.f
  abs.hm.ug        <- exp.hmouth.tot * ingest.abs.f
  abs.inhal.ug     <- cb$dose.inhal.tot * cb$inhal.abs.f
  abs.tot.ug       <- abs.dermal.ug + abs.ingest.ug + abs.inhal.ug
  abs.tot.mgkg     <- abs.tot.ug/(1000*cb$weight)
  urine.tot.ug     <- abs.tot.ug * cb$urine.f

  fexp <- as.data.table(cbind(cb, rem.bath, rem.wash, rem.hmouth,
                              rem.brush, rem.absorb, exp.hmouth.tot, dose.intake.ug,
                              dose.intake.mgkg, abs.dermal.ug, abs.ingest.ug,
                              abs.inhal.ug, abs.tot.ug, abs.tot.mgkg, urine.tot.ug,abs.hm.ug))
  return(fexp)
}


# Notes for SHEDS-HT main module
#
# Notes for run()
#   Arguments: run.file, wd
#              The user specifies the name of the run.file for this run.
#              Many different run.files may be set up for special purposes.
#              The one being invoked must be present in the /inputs folder.
#              No keyword is needed, but the file name must be quoted.  No
#              file extension is needed, and .CSV is assumed.
#              The wd value should be set once for each SHEDS-HT installation
#              by replacing the default value in the definition of run().
#
#  setup()  This function loads the required R packages and sources the
#           modules.  The user might need to download the packages if they
#           are not already present.
#
#  specs    read.run.file() reads the specified run.file, which must be a .CSV
#           file in the /inputs folder.  This file contains one header line,
#           followed by lines each containing one keyword and a value,
#           separated by an equal sign.  A default file called 'run.file.csv'
#           contains an example for each valid keyword.  The 'specs' object
#           is the only global variable in SHEDS-HT. It is not altered after
#           it is defined, and can be accessed by any function in SHEDS-HT.
#
# act.diaries There are two types of diaries in SHEDS-HT.  The first are the
#             activity diaries, which indicate the amount of time and level
#             of metabolic activity in various 'micros'.  Each line of data
#             represents one person-day (24 hours).
#
# chem.props  These are the chemical properties required for SHEDS-HT. The
#             default file was prepared from publically available databases
#             using a custom program (not part of SHEDS-HT).  The default file
#             contains about 17 numerical inputs per chemical, but most are
#             not used. The required propoerties are molecular weight (MW),
#             vapor pressure (vapor), solubility (solub), octanol-water
#             partition coefficient (KOW), air decay rate (decay.a), decay
#             rate on surfaces (decay.s), and permeability coefficient (Kp).
#             Additional variables may be present and are ignored.
#
# master.chem A condensed version of chem.props that contains only the
#             modified CAS numbers (with underscores) and a frequency count
#             for each.  All frequency counts will be one, unless the same
#             CAS number appears multiple times in the database (which could
#             signify a problem).  This list of CAS numbers is useful for
#             searching, when the entire chem.props object is cumbersome.
#
# diet.diaries  These are daily diaries of dietary consumption by food group.
#               Each line represents one person-day, with demographic
#               variables followed by amounts (in grams) for a list of food
#               types indicated by a short abbreviation on the header line.
#               The user may supply their own food consumption database with
#               a customized list of food types, however, the same food types
#               need to occur on the Scenarios file for any exposure to result.
#
#  foods  This is a list of the labels for the food categories on the
#         consumption database.  An exposure scenario can generate exposure
#         only if its food type is in this list.
#
# exp.factors   This data set contains the distributional parameters for
#               the exposure factors, which include the fraction of chemical
#               available for transfer during dermal contact (avail.f), the
#               probability of taking a bath or shower on a given day (bath.p),
#               the dermal transfer coefficient (dermal.tc), the fraction of
#               hand area (noth hands) typically mouthed (handmouth.area.f),
#               hound mouthing frequency per hour (handmouth.freq), hand
#               washing frequency per day (handwash.freq), the absorption
#               fraction for inhaled chemcial (inhal.abs.f), the exposure
#               ratio for object-to-mouth contact (om.ratio), the fraction
#               of dermal loading removed by a bath or shower (rem.bath.f),
#               the fraction of dermal loading removed by 'other' processes,
#               including brushing, rubbing, falling off, touching objects
#               and surfaces, and shedding of skin (rem.brushoff.f), the
#               fraction of dermal loading absorbed through the skin
#               (rem.dermal.abs.f), the fraction of dermal loading removed
#               by hand washing (rem.handwash.f), and the fraction of
#               absorbed dose excreted in urine (urine.f).
#               All of these variables may have age or gender-dependent
#               distributions, although in the absence of data, many are
#               assigned a single distribution from which all persons are
#               sampled.
#
# fug.vars  The fugacity input file contains distributions for 36 variables
#           needed to support the in-home fugacity modeling.  Two scenario
#           types need these variables: indir.ug models a single large
#           application of chemical, while indir.y0 models a continual
#           slow release from a chemical-containing object.
#
# media     This is a list of the potential contact media for the model.
#           Each is found in a specific microenvironment (micro). Contact.p
#           gives the probability of contact (that is, exposure) with the given
#           media type during the time spent in that micro.
#
# media.sur   A list of surface media
#
# media.air   A list of air media
#
# physio    The three physiological variables of interest are weight, height,
#           and body mass index.  The input file contains regression
#           parameters for each of them, for various age and gender groups.
#
# pop       The population file contains counts by gender and each year of
#           age from the 2000 U.S. census.  When a large age range is
#           modeled, this ensures that SHEDS chooses age and gender with the
#           correct overall probability.
#
# all.scenarios This is the master list of all chemicals, and all exposure
#               scenarios specific to each chemical, to be evaluated in the
#               current model run.  The user may create multiple scenarios
#               files for special purposes, for example, for selected chemical
#               classes, or selected exposure pathways. A model run consists
#               of two nested loops: the outer loop over chemicals and the
#               inner loop over the scenarios specific to that chemical.
#
# act.pools The activity diaries are assigned to "pools" which reflect the
#           age, gender, weekend, and season assignments on the diaries. The
#           last three variables must match the simulated day, and age must
#           agree within a percentage given by specs$age.match.pct.  For
#           example, a 35 year old person with age.match.pct=20 may be assigned
#           a diary with a nominal age between 28 and 42, inclusive. Act.pools
#           contains a list of suitable diaries for each year of age, for
#           each gender, weekend, and season combination.
#
# diet.pools  The dietary consumption diaires are pooled similarly to the
#             activity diaries, but there are fewer pools because only age
#             and gender are used (not weekends or seasons).
#
# gen.facs    These are the general exposure factors, meaning those that are
#             are not media-specific. However, these exposure factors are
#             potentially gender and/or season specific.  Each possibility is
#             assigned a row number pointing to the approproate distribution
#             on the exp.factors data set.  Gen.facs has 8 rows per variable
#             (2 genders x 4 seasons), regardless of the number of different
#             distributions used for that variable.
#
# med.facs    These are the media-specific exposure factors.  The present
#             version of the model has only 3, but the code will accept more.
#             The default media.csv file flags only 3 media as being potential
#             exposure locations, so each variable has 24 rows (2 genders,
#             4 seasons, 3 media) when all ages share the same distribution.
#             If a variable has N age categories (each with its own
#             distribution) then there are (24 N) rows for that variable.
#
# chem.scen   This data set has a list of chemicals, and the scenarios for
#             that chemical that have potential for exposure.  It is possible
#             for one chemical to have two or more "dermal" scenarios,
#             provived each is in a different "category".  Each combination
#             of chemical, category, and scenario must be unique. The last
#             cariable on chem.scen is "count", which is the number of rows
#             of data on all.scenarios devoted to this combination.
#
# nscen       Nscen lists the number of scenarios for each chemical. The
#             number of rows on nscen equals the number of chemicals.
#
# last.chem   The run stops when this chemical number has been finished.
#             It is the lesser of the number of rows on nscen and the
#             specs$maxChem setting.  Note that if specs$firstChem is not 1,
#             then the number of chemicals actually processed is less than
#             last.chem.
#
# base    The base data set contains all the chemical-independent information
#         needed for the exposure assessment.  There is one row for each
#         simulated person.  The variable include age, gender, weight, diet
#         and activity diaries, food consumption, minutes in each micro, and
#         the evaluation of the exposure factors for that person.  This
#         data remains constant across all chemicals and scenarios.  Once the
#         base data set is ready, a pair of nested loops over chemicals and
#         scenarios are entered.
#
# fexp    The first step inside the loop over chemicals is to remove the
#         fexp data set, if one exists.  Fexp stands for "final exposure"
#         and it contains all the information on exposure to one chemical.
#         When starting a new chemical, the old one is erased to ensure that
#         results from different chemicals cannot be accidentally combined.
#
# chem    This is the modified CAS number of the current chemical.  In SHEDS
#         the cas number contains underscores "_" rather than dashes "-".
#
# cdata   The list of scenario-specific information for the current chemical.
#         This is much more detailed than chem.scen or cs.
#
# cs      The scenarios specific to this chemical.
#
# ns      The number of scenarios on cs.
#
# crow    The row number on the chem.props and master.chem objects. If crow
#         is missing, then a message is printed and this chemical is not
#         considered further.
#
# cb      The cb object is a copy of "base" with columns added for exposure
#         variables.  These are initialized to zero, in case later operations
#         do not change them.  For summarization purposes, it is better to
#         have zeroes than to have undefined quantities.
#
# chemical  This is the name of the chemical (as opposed to CAS number).
#
# ns.good This is a running count of the number of valid scenarios for this
#         chemical.  The chemical will be summarized later if ns.good>0.
#
# scen    This is a codeword for the type of scenario. It appears on each
#         row of the scenarios input file. The acceptable scenario types are
#         seen in a series of IF statements following shortly in the code.
#
# categ   The catagory for this scenario. It has little function in SHEDS,
#         except to allow multiple scenarios of the same type for the same
#         chemical.  For example, direct dermal exposure "dirderm" is a
#         scenario type which might have multiple instances for the same
#         chemical, perhaps in pesticide use, paint, cleaners, and so on.
#
# bad     A flag indicating the status of the variables for this scenario.
#         If bad=0 then all required variables have been defined.  Otherwise,
#         the scenario is skipped and the value of bad may indicate what
#         was wrong.
#
# sdata   The chemical-scenario data specific to this combination of chemical
#         and scenario.  For the DIETARY scenario, a call to checkfoods is
#         made to ensure that the same food groups do not appear more than
#         once, which can sometimes happen and would cause code failure.
#
# factors These are the distributional parameters for the variables in sdata.
#         This is where the instructions from the scenarios input file are
#         parsed and interpreted.  For example, here the prevalence is
#         converted from a percentage to a fraction. Also, the parameters
#         for lognormal distributions are converted to SHEDS parameters.
#         These operations are performed once for each distribution, not
#         repeatedly for each person.
#
# indices A mapping showing the appropriate row of the "factors" object for
#         each combination of variable, gender, and age.
#
# sf      These are the person-specific values for the randomly sampled
#         variables from sdata. As with other randomly sampled variables
#         in SHEDS, all samples are independent.  Each row of sf contains
#         the values for one person.  Following the evaluation of sf, the
#         code contains a series of IF statements that process the various
#         exposure scenarios.  Each scenario will match one of these.
#
# dirderm   The dir.dermal function calculates direct dermal exposure that
#           occurs when a chemical-containing product is used.  This does
#           not include later contact with treated objects, which is
#           indirect exposure.
#
# diringest The dir.ingest function considers chemicals ingested during or
#           immeidately after use. For example, chemcials found in
#           toothpaste.  However, food and drinking water are not included.
#
# dirinhaer This covers inhalation from aerosols, like hairspray and similar
#           products that are directly injected into the air on or around
#           the exposed person.
#
# dirinhvap This covers inhalation of vapor from volatile chemicals.
#
# dietary   This covers chemicals found in food and drinking water.
#
# downthedrain  This scenario is different from the others because the
#               target of the exposure is not the simulated person, but
#               the sewer system. This type of scenario gained widespread
#               notice in the 1970s with the issue of phosphates in laundry
#               detergent.
#
# indir.ug   There are two indirect scenarios in SHEDS at present. This
#             one models the spread of chemical around a house after a
#             single application or treatment.  The time since the
#             treatment occurred is an important variable.
#
# indir.Y0    This is the second of the indirect scenarios in SHEDS. It
#             examines the effect of constant emission sources.  For
#             example, carpets and furniture may give off small amounts of
#             formaldyhe indefinitely.
#
# fexp      After the specfic exposure scenarios are evaluted, SHEDS calls
#           the post.exposure function.  The results are written to fexp.
#           The post.exposure function considers the fate of dermal loading,
#           by evaluating washing, bathing, hand-to-mouth transfer, dermal
#           absorption, and brush-off terms.  It also calculates absorption
#           for all pathways and the amount of chemical in the urine.
#
# summarize.chemical   If a chemical has scenarios that have been evaluated,
#                     the summarize.chemical function is called to write
#                     an output file specific to that chemical. This file
#                     is named "CAS_....", where the dots are the CAS number.
#                     The file is in the folder indicated by the runName.
#
# setup()   This function is called at the start of a run. It sources all
#           the SHEDS code and loads the required R packages.  It also
#           sets the working directory for input and output files.
#
# act.diary.pools() The act.diary.pools function assigns the activity diaries
#                 to pools. A given diary may belong to many pools, since
#                 every year of age has its own pool.  The output of this
#                 function is "pool", which is vector of lists.  The length
#                 of "pool" is the product of #genders, #seasons, and #ages.
#                 Each element is a list of acceptable diary numbers.
#                 For example, pool[100] may have the name M0P99, which
#                 indicates that it is for males, on weekdays, in spring,
#                 for age=99.  Large pools may contain a list of several
#                 hundred diary numbers.  The code contains four loops and
#                 on each step perform a subsetting of the list of diary
#                 numbers.
#
# diet.diary.pools()  This is like act.pools, but for dietary consumption.
#                   The pool types are referenced by gender and age only,
#                   so there are usually fewer than for act.pools.
#
# gen.factor.tables() This function constructs the tables of general exposure
#                   factors, meaning the ones that are not media-specific.
#                   The input comes from the exposure factors file. The
#                   code loops over the input rows, making multiple copies
#                   of each with the approprate age, gender and season. For
#                   simplicity these new rows are tacked on to a copy of
#                   the original data set. At the end, the original rows
#                   are deleted.
#
# med.factor.tables() This is similar to gen.factor.tables, except that it
#                   operates on the three media-specific variables in
#                   exp.factors, namely avail.f, dermal.tc, and om.ratio.
#                   This function contains an extra nested loop over media
#                   when compared to gen.factor.tables.
#
# chem.scenarios()   This function summarizes the all.scenarios dataset by
#                   condensing each chemical-scenario combination into one
#                   row. The number of rows of data for this combination on
#                   all.scenarios is recorded here.
#
# select.people()  This is the first real step in the modeling process. It
#                 first fills an array "q" with uniform random numbers,
#                 with ten columns because there are 10 random variables
#                 defined by this function. There are "n" rows, which is
#                 the number of persons, capped at specs$set.size (typically
#                 5000).  Gender is selected from a discrete (binomial)
#                 distribution where the counts of males and females in the
#                 study age range determines the gender probabilities. Age
#                 next, separately for each gender. The counts by year of
#                 age is chosen for the appropriate gender as used as
#                 selection weights. Season is assigned randomly (equal
#                 weights) using the ones on specs$seasons. "Weekend" is
#                 set to one or zero, with a chance of 2/7 for the former.
#                 All the above are vector operations, with an independent
#                 outcome based on "q" for each simulated person.
#
#                 The next block of code assigns physiological variables.
#                 "r" is a vector of offsets in various age-gender lists,
#                 with one value for each person.  Weight is lognormal in
#                 SHEDS, so a normal is sampled first and then exp() is
#                 applied.  This means that the weight parameters refer to
#                 the properties of log(weight), which were fit by linear
#                 regression. For reasons I do not recall, height is
#                 calculated in two ways, and the larger value is used.
#                 Bmr is basal metabolic rate, which uses a regression.
#                 A minimum bmr is set to prevent extreme cases from
#                 becoming zero or negative.  Bva is the alveolar breathing
#                 ventilation rate corresponding to bmr. The SHEDS logic
#                 sets activities in each micro to be a multiple of these
#                 rates, with outdoor rates higher than indoor ones, which
#                 are higher than sleep. This affects the inhaled dose.
#                 Surfarea is the skin surface area, which uses regrssions
#                 based on height and weight for 3 age ranges.  To avoid
#                 using IF statements or FOR loops over thousands of people,
#                 all 3 are calculated for everybody, and the appropriate
#                 one is chosen by subsetting.
#
#                 The next step is to assign diaries.  Here, a FOR loop over
#                 persons is used.  The values "atype" and "dtype" are the
#                 names of the diary pools appropriate for each person. The
#                 call to "Distrib" picks one from an empirical distribution
#                 consisting of the list of diary numbers for that pool. The
#                 "people" data set consists of all the assigned values. The
#                 final step is to retrieve the actual data from the chosen
#                 activity and diet diaries, and the result becomes "pd".
#
# add.media()  First, an array "q" of uniform random samples is generated
#             with rows=#people (in this set) and cols=#potential exposure
#             media. The array "mic" is the number of minutes on the activity
#             diary in the relevant micro, with a row for each person and a
#             column for each of the exposure media. The array "dur" is
#             similar except that each duration is multiplied by the
#             relevant "contact.p".
#
# add.factors()  This function adds the specific exposure factors to "pdm".
#               Here, "specific" means taking into account the age, gender,
#               season, and sometimes the exposure medium, for each person.
#               First, "w" is determined, which is the number of general
#               factors plus the product of the number of media-specific
#               factors and the number of surface media.  Air media do not
#               have media specific factors in this version of SHEDS. An
#               array "q" of uniform random samples is generated, one row
#               per person and "w" columns.  A zero matrix "r" is defined,
#               of the same size.
#
#               The media-specific factors are handled first. Two nested
#               loops over variable and surface type generate the values,
#               which are stored in "r".  Then another loop covers the
#               general factors. The "p" data set contains the age, gender,
#               and season for each person, and the "join by" operation
#               combines them. The actual evaluation is handled by
#               eval.factors(), discussed below.
#
#               One of the exposure factors is handwash.freq. This was also
#               part of SHEDS-Multimedia, where technically it represented
#               the mean number of hours in the day with hand washing events.
#               An important aspect of that model was that each person was
#               followed longitudinally, and the actual number of hand washes
#               on each day varied from one day to the next. Because of this,
#               the distribution for handwash.freq did not need to be
#               restricted to integer values, as (for example) a mean of 4.5
#               per day is acceptable and achievable, while choosing integer
#               numbets of hand washes each day.  One of the early goals
#               with SHEDS-HT was to attempt to reproduce selected results
#               from SHEDS-Multimedia. Therefore, similar logic was built
#               into the current model. The hand.washes variable is sampled
#               from a distribution centered on handwash.freq, and then
#               rounded to the nearest integer.
#
#               The "bath" variable is another difficult concept. In theory,
#               baths and showers are recorded on the activity diaries. In
#               practice, the CHAD activity database was constructed from
#               around 20 separate studies, and some of those did not contain
#               enough detail to identify separate bath or shower events. The
#               result is that about half of all diaries record such events,
#               but the true rate in the population is higher.  The bath.p
#               variable was created to address this. It represents the
#               proability that a non bath/shower activity diary should
#               actually have one. Therefore, if the diary has one, then
#               SHEDS automatically has one. Otherwise, a binomial sample
#               using bath.p as the probability is drawn. A bath/shower occurs
#               unless both of these are zero.
#
#               The effectiveness of hand washes or bath/shower at removing
#               chemical from the skin is determined in post.exposure().
#
# eval.factors() This function is rather abstract, so it may be hard to follow.
#               The three arguments to eval.factors() are the vector of row
#               numbers (one per person), the vector of random samples, and
#               the list of distributions for the exposure factors. A list
#               "ur" of the unique rows that are actually used is created, to
#               avoid evaluating unused distributions. The FOR loop covers
#               each of these distributions separately. The index "j" is the
#               actual row number on exp.factors (here called "ef"). For
#               example, suppose for a given variable all the values of p$row
#               (here called "r") were 4,5, or 6.  Then the FOR loop is run
#               through 3 times.  All the rows on exp.factors other than
#               4, 5, or 6 are skipped, usually because they apply to other
#               variables.  For j=4, the key line in the code populates all
#               rows of the output vector z for which r=4 with samples from
#               the "distrib" function called with the parameters from row
#               "j" of "ef".  The number of samples returned by distrib
#               is given by the length of the "q" arguments, which is the
#               list of desired quantiles.  Setting this to q[r==j] subsets
#               the vector of quantiles (originally one per person) to a
#               vector of just those for people assigned to row "j" for this
#               exposure factor.  In this case, "distrib" does not generate
#               any randomness, it just applies deterministic transforms to
#               the "q" values.  This is a slow step in SHEDS, and by doing
#               only the required persons and required rows of exp.factors,
#               the program runs faster. An explicit loop over persons and
#               variables would be clearer, but slower, because each call to
#               distrib would return just one value.  For a large run, the
#               current method may return thousands of values per call, which
#               is much more efficient.
#
# make.cbase    This function extends the "base" data set to include a set
#               of exposure variables specific to the current chemical. At
#               this point, all these new variables are initialized to zero.
#
# create.scen.factors() This function takes the information on distributions
#                     from the all.scenarios data set (which comes from the
#                     "scenarios" input file) and converts it into the
#                     parameter set needed by SHEDS. At present, these steps
#                     are 1) converting prevalence from a percentage to a
#                     binomial form, 2) converting CV to standard deviation
#                     for normals, and 3) for lognormals, converting mean and
#                     CV to par1 and par2.  Also, some variables are renamed.
#
#                     Note that "prevalence" in SHEDS becomes a binomial
#                     distribution, which returns a value of either 0 or 1
#                     when evaluated. Each simulated person either "does" or
#                     "does not" partake in this scenario. Similar logic
#                     applies to the "frequency" variable, except that the
#                     returned values may be larger than one (that is, 2 or
#                     more) for very frequent scenarios.  All the exposure
#                     equations contain the "prevalence" variable.  If is set
#                     to one for that person, the exposure is as expected,
#                     but if prevalence=0 then no exposure occurs.
#
# scen.factor.indices() This function resembles gen.factor.tables(), except that
#                     it applies to scenario-specific factors that come from
#                     the scenarios input file. The difference is that these
#                     factors do not allow seasonal variation, so one of the
#                     loops is missing. Also, these factors are not media-
#                     specific, altough they obviously are scenario-specific,
#                     and most scenarios include just one surface medium.
#
#                     This function, and eval.sscen.factors() below, are
#                     evaluated separately from the other factors because
#                     these may change with every chemical and scenario (so
#                     they are constantly being updated), whereas the other
#                     factors are evaluated once overall for each person.
#
#
# dir.dietary()   This is the first of the specific exposure scenarios, for
#                 the "dietary" pathway, which means chemicals found in foods
#                 (including water and drinks). The data exist in two parts:
#                 one is the food consumption database which has intake in
#                 grams per day, and the other is the scenarios file which
#                 lists the chemical content in certain food groups. It is
#                 critical that the food groups on these two files are
#                 consistent. They do not have to have full overlap because
#                 any food group not mentioned on the scenarios file is
#                 assumed to be free of the chemical in question. Also, if
#                 the chemical is present in an unconsumed food group, no
#                 exposure results.  However, all food groups in common to
#                 both files should be defined similarly.
#
#                 In this function, "foods" is a list of the names of the
#                 food groups, often two letters in length, but may be up to
#                 seven letters. The FOR loop picks the name of each food
#                 group as a string, and converts it to a variable name. For
#                 example, the food group "FV" may have corresponding
#                 variables residue.FV, zeros.FV, and nonzeros.FV on the
#                 cbase data set. If these variables are not present, no
#                 exposure results.  The variables may be present, but an
#                 individual person may still receive zero exposure because
#                 not all samples of that food are contaminated, as indicated
#                 by the "zeros" value.  THe "nonzeros" value is used to
#                 determine the likelihood of contamination, with the residue
#                 variable determining the amound found (when it is nonzero).
#                 The dietary exposure is the product of the consumption
#                 (in grams) and the residue (in micrograms of chemical per
#                 gram of food), summed over all food groups.  Three other
#                 vectors are returned (dermal exposure, inhalation exposure,
#                 and inhalation dose), but these are currently set to zero
#                 for the dietary pathway. In principle, eating food could
#                 result in dermal exposure (for finger foods), or inhalation
#                 (as foods may have noticeable odors), but these are not
#                 large enough to be of concern at present.
#
# dir.dermal()     The dermal exposure scenario is relatively straightforward.
#                 The prevalence reflects the fraction of the population who
#                 use this scenario at all. The frequency is the mean number
#                 of times per year this occurs among the "doers".  Since
#                 SHEDS is considering one random day, the frequency is
#                 divided by 365 and then passed to the p.round (probabilistic
#                 rounding) function, which rounds either up or down to the
#                 nearest integer. Very common events may happen more than
#                 once in a day. The "mass" variable refers to the  mass of
#                 the product in grams in a typical usage event. The
#                 composition is the percentage of that mass that is the
#                 chemical in question.  The final variable is "resid", which
#                 measures the fraction that is likely to remain on the skin
#                 when the usage event ends.  The dermal exposure is the
#                 product of the above variables, with units converted to
#                 micrograms by the final factor of 1E6.
#
# dir.ingested()   This scenario is for accidental ingestion during product
#                 usage. Typical examples are toothpaste, mouthwash, lipstick
#                 or chap stick, and similar products used on the face or
#                 mouth. The variables are generally similar to the dermal
#                 scenario, including prevalence, frequency, mass, and
#                 composition. The distinct variable is "ingested", which is
#                 the precentage of the mass applied that becomes ingested.
#                 Since these products are not intended to be swallowed, this
#                 should typically be quite small (under 5%). The final factor
#                 in the exposure equation is a units conversion from grams
#                 to micrograms.
#
# dir.inh.aer()     This scenario considers inhalation exposure from the use of
#                 aerosol products. Typicaly examples include hairspray and
#                 spray-on mosquito repellent. The prevalence, frequency,
#                 mass, and composition variables have the same meaning as
#                 for the dir.dermal and dir.ingested scenarios.  Unlike those
#                 scenarios, "exposure" for the inhalation pathway has units
#                 of micrograms per cubic meter, reflecting the average air
#                 concentration of the chemical. An "airconc" variable is
#                 defined using mass, composition, frac.aer, and volume.
#                 Frac.aer is the fraction of the product mass that becomes
#                 aerosolized, and the volume affected by the use is
#                 approximated to allow the calculation of a concentration
#                 or density. Defaults are set in the code if these variables
#                 are missing from the input file. Since "exposure" depends
#                 on the time-averaged concentration, a duration is necessary.
#                 For example, if one spends five minutes in an aerosol cloud
#                 and the rest of the day in clean air, the daily exposure is
#                 the cloud concentration multiplied by 5/1440 (where 1440 is
#                 the number of minutes in a day).
#
#                 dirinhAer() also calculates the inhaled dose, in units of
#                 micrograms per day. The dose equals the product of exposure
#                 (g/m3), basal ventilation rate (m3/day), the METS factor
#                 (typically people inhale air at an average of 1.75 times
#                 the basal rate to support common daily activities), and a
#                 conversion factor of 1E6 from grams to micrograms.
#
# dir.inh.vap()     This is another inhalation scenario, but this one is for
#                 vapor. For example, painting will result in the inhalation
#                 of vapor, but it does not involve aerosols (unless it is
#                 spray paint).  For this scenario, the vapor pressure and
#                 the molecular weight are relevant variables for determining
#                 exposure.  These do not have to be entered on the scenarios
#                 file, they are looked up from the chem.props file. That file
#                 reports vapor pressure in Pascals; SHEDS converts it to
#                 atmospheres by dividing by 1E5, and then caps it at one. As
#                 a practical matter, the vapor pressure in an open area will
#                 not exceed one atmosphere.
#
#                 As indicated in the function itself, "evap" is effectively
#                 an evaporated mass, calculated using the product mass,
#                 composition (converted from percent to a fraction), the
#                 vapor pressure as a surrogate for partial pressure, and
#                 duration of product use.  The duration term is made unitless
#                 by dividing by 5 (minutes), which is an assumed time
#                 constant. The logic is ???
#
#                 The maximum allowed air concentration is "maxconc", which
#                 is the point at which evaporation ceases (or technically,
#                 is balanced by condensation).  For chemicals used for a
#                 short duration, or with low vapor presssure, this point
#                 might never be reached before usage stops.
#
#                 Once the effective air concentration is established, the
#                 rest of theis function handles exposure just like the
#                 dir.inh.aer() function.
#
# down.the.drain.mass()  This may be the simplest of all the current scenarios.
#                     It evaluates the amount of chemical entering the waste
#                     water system, on a per person-day basis.  The "exposure"
#                     is to a system, not a person, but this method uses one
#                     person's actions to estimate their contribution to the
#                     total. The amount is the product of the prevalence, the
#                     frequency (converted to an integer count using the
#                     p.round function), the product mass, the composition
#                     (converted from a percentage to a fraction), and the
#                     fraction going down the drain (also converted from a
#                     percentage).  The result is in grams per person-day.
#
# indir.exposure()   Indirect exposure happens after a product is no longer
#                   being used or applied, due to chemical lingering on
#                   various surfaces or in the air. People who come along
#                   later may get dermal or inhalation exposure from residual
#                   chemical in the environment.
#
#                   SHEDS-HT currently has two indirect exposure scenarios.
#                   "INDIR.FUG" applies to a one-time chemical bolus applied
#                   to a house, whereas "INDIR.Y0" applies to continual
#                   releases. Both scenarios consist of two parts: the first
#                   determines the appropriate air and surface concentrations,
#                   That code is in the Fugacity module.  The second part is
#                   the exposure calculation by indor.exposure(). Both types
#                   of indirect exposure scenarios call this function.
#
#                   The surface and air concentrations from the fugacity code
#                   are premised on the product use actually occurring. Hence,
#                   indir.exposure() starts by multiplying those concentrations
#                   by the prevalence (which is either 0 or 1, evaluated
#                   separately for each person).  The affected media are
#                   assumed to be "inawkair" and "inawksur", so the contact
#                   durations for these media are used in the calculations.
#                   For air, the exposure is the average daily concentration,
#                   which is the event concentration multiplied by the fraction
#                   of the day spent in that event. The inhaled dose is the
#                   product of the exposure, the basal ventilation rate, and
#                   the PAI factor (multiplier for the basal rate). A factor
#                   of 1E6 converts the result from grams per day to micrograms
#                   per day.
#
#                   Dermal exposure results from skin contact with surfaces.
#                   The surface concentration (ug/cm2) is multiplied by the
#                   fraction available for transfer (avail.f, unitless), the
#                   transfer coefficient (dermal.tc, in cm2/hr), and the
#                   contact duration (hr/day).  The result is the amount of
#                   chemical transferred onto the skin (ug/day).
#
# post.exposure()  This function resolves the fate of chemical after initial
#                 exposure has occurred.  The same function applies to all
#                 exposure scenarios. The dermal exposure is the most
#                 complicated, as there are five removal methods. All five
#                 are randomly sampled.  While the sum of the five means is
#                 close to one, the sum of five random samples might not be,
#                 so these samples are treated as fractions of their sum.
#
#                 The cb$bath variable is either 0 or 1, so the bath removal
#                 term is either zero or the sampled bath removal efficiency.
#                 The brushoff term is also simple. The handwashing and hand-
#                 to-mouth transfer terms are non-linear, because higher
#                 frequencies have less chemical available for removal on
#                 each repetition.  The functions in place were fitted to
#                 output from SHEDS-Multimedia, which were summed to daily
#                 totals.
#
#                 The final removal term is dermal absorption. The base value
#                 is multiplied by the Kp factor from the chemical properties
#                 database, and divided by the value for permethrin, as the
#                 values from SHEDS-Multimedia to scale the other removal
#                 terms were based on a permetrin run.
#
#                 The five terms are evaluated separately for each person, as
#                 is their sum. Each is then converted to a fraction of the
#                 whole. The fractions may be quite different from one person
#                 to another.  FOr one, perhaps 80% of the dermal loading is
#                 removed by a bath/shower, while for another person it is
#                 0% because they did not take one.  For the latter person,
#                 other four removal terms are (on average) five times larger
#                 than for the former person, because together they account
#                 for 100% of the removal, instead of just 20%.
#
#                 The rest of the post.exposure function is mostly a matter of
#                 bookkeeping. The hand-to-mouth dermal removal term becomes
#                 an ingestion exposure term. Summing exposures across routes
#                 is dubious, in part because inhalation exposures use
#                 different units from the others, and because much of the
#                 dermal exposure never enters the body, and finally because
#                 summing dermal and ingestion exposures may double-count
#                 the hand-to-mouth term. However, intake dose may be summed.
#                 In SHEDS, "intake dose" is the sum of the inhaled dose
#                 (which is the amount of chemical entering the lungs in
#                 ug/day), the ingestion exposure (which is the amount
#                 entering the GI tract in ug/day), and the dermal absorption
#                 (the amount penetrating into or through the skin, so it
#                 cannot otherwise be removed, in ug/day). The "absorbed dose"
#                 is also calculated: for dermal it is the same as the intake
#                 dose, but for ingestion and inhalation there is another
#                 absorption factor, which was set on the exposure factors
#                 input file. An estimate of the chemical in urine (in ug/day)
#                 is made.  Both the intake dose and the absorbed dose are
#                 reported in both (ug/day) and in (mg/kg/day). Note that the
#                 latter requires the body weights of each individual. These
#                 cannot be obtained from the former just by knowing the
#                 average body weight in SHEDS.
#
#                 The above variables for each simulated person are written
#                 to the "fexp" object (the name stands for "final exposure").
#                 Each chemical writes over the previous fexp, so the data
#                 must first be summarized and written to an output file.
