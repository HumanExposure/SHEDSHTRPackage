#' read.run.file
#'
#' Each SHEDS run has its own "run.file" that the user prepares before the run.   This file contains all the settings and file references needed for the run.  Read.run.file occurs at the start of each SHEDS run.  The contents of the run.file are examined, and stored in the R object "specs".
#' 
#' @export

read.run.file = function(run.file="run_test.txt") {

  if(exists("specs"))   rm(specs,inherits = TRUE)
  run.name         <- "prods"
  n.persons        <- 10
  person.output    <- "yes"
  source.output    <- "yes"
  min.age          <- 3
  max.age          <- 5
  genders          <- c("M","F")
  seasons          <- c("P","S","F","W")
  details          <- 1
  age.match.pct    <- 20
  run.seed         <- 876144637
  set.size         <- 1000
  act.diary.file   <- "Activity_diaries.csv"
  chem.props.file  <- "Chem_props.csv"
  diet.diary.file  <- "Diet_diaries_0.csv"
  exp.factor.file  <- "Exp_factors.csv"
  fugacity.file    <- "Fugacity.csv"
  media.file       <- "Media.csv"
  physiology.file  <- "Physiology.csv"
  population.file  <- "Population.csv"
  source.vars.file <- "Source_vars.csv"
  source.scen.file <- "Source_scen_prods.csv"
  source.chem.file <- "Source_chem_prods.csv"
  chem.list        <- ""
  n.chem           <- 0

  if (!grepl("[.]txt$", run.file) & !grepl("[.]csv$", run.file)) {
     run.file = paste0(run.file,".txt") }
  a <- read.table(paste0("inputs/",run.file),skip=1,sep="=",
                  stringsAsFactors=FALSE,strip.white=TRUE)
  a$V1 <- tolower(a$V1)
  a$V2 <- gsub(" ","",a$V2)


  for (i in 1:nrow(a)) {
    if (str_trim(a$V1[i])=="run.name")         run.name         <- str_trim(a$V2[i])
    if (str_trim(a$V1[i])=="n.persons")        n.persons        <- as.integer(a$V2[i])
    if (str_trim(a$V1[i])=="person.output")    person.output    <- str_trim(a$V2[i])
    if (str_trim(a$V1[i])=="source.output")    source.output    <- str_trim(a$V2[i])
    if (str_trim(a$V1[i])=="min.age")          min.age          <- as.integer(a$V2[i])
    if (str_trim(a$V1[i])=="max.age")          max.age          <- as.integer(a$V2[i])
    if (str_trim(a$V1[i])=="genders")          genders          <- c(str_trim(a$V2[i]))
    if (str_trim(a$V1[i])=="seasons")          seasons          <- c(str_trim(a$V2[i]))
    if (str_trim(a$V1[i])=="details")          details          <- as.integer(a$V2[i])
    if (str_trim(a$V1[i])=="age.match.pct")    age.match.pct    <- as.numeric(a$V2[i])
    if (str_trim(a$V1[i])=="run.seed")         run.seed         <- as.integer(a$V2[i])
    if (str_trim(a$V1[i])=="set.size")         set.size         <- as.integer(a$V2[i])
    if (str_trim(a$V1[i])=="act.diary.file")   act.diary.file   <- str_trim(a$V2[i])
    if (str_trim(a$V1[i])=="chem.props.file")  chem.props.file  <- str_trim(a$V2[i])
    if (str_trim(a$V1[i])=="diet.diary.file")  diet.diary.file  <- str_trim(a$V2[i])
    if (str_trim(a$V1[i])=="exp.factor.file")  exp.factor.file  <- str_trim(a$V2[i])
    if (str_trim(a$V1[i])=="fugacity.file")    fugacity.file    <- str_trim(a$V2[i])
    if (str_trim(a$V1[i])=="media.file")       media.file       <- str_trim(a$V2[i])
    if (str_trim(a$V1[i])=="physiology.file")  physiology.file  <- str_trim(a$V2[i])
    if (str_trim(a$V1[i])=="population.file")  population.file  <- str_trim(a$V2[i])
    if (str_trim(a$V1[i])=="source.vars.file") source.vars.file <- str_trim(a$V2[i])
    if (str_trim(a$V1[i])=="source.scen.file") source.scen.file <- str_trim(a$V2[i])
    if (str_trim(a$V1[i])=="source.chem.file") source.chem.file <- str_trim(a$V2[i])
    if (str_trim(a$V1[i])=="chemical") {
        n.chem <- n.chem+1; chem.list[n.chem] <- str_trim(a$V2[i]) }
  }

  dir  <- paste0("output/",run.name)
  if(!file.exists(dir)) dir.create(dir,recursive=TRUE)
  set.seed(run.seed)
  p.out <- tolower(substr(person.output,1,1))
  person.output <- ifelse (p.out=="y" | p.out=="1", 1, 0)
  s.out <- tolower(substr(source.output,1,1))
  source.output <- ifelse (s.out=="y" | s.out=="1", 1, 0)
  if(n.chem==0) chem.list <- NULL else {
    chem.list <- sort(unique(trimzero(str_replace_all(chem.list,"-","_"))))
    n.chem <- length(chem.list)
  }
  gend <- toupper(genders)
  b    <- ""
  if (str_detect(gend,"M"))  b <- c(b,"M")
  if (str_detect(gend,"F"))  b <- c(b,"F")
  genders <- b[2:length(b)]
  seas    <- toupper(seasons)
  b       <- ""
  if (str_detect(seas,"P"))  b <- c(b,"P")
  if (str_detect(seas,"S"))  b <- c(b,"S")
  if (str_detect(seas,"F"))  b <- c(b,"F")
  if (str_detect(seas,"W"))  b <- c(b,"W")
  seasons <- b[2:length(b)]
  if (n.persons<=0)        stop("\n No persons to model \n")
  if (is.null(genders))   stop("\n No genders selected \n")
  if (is.null(seasons))   stop("\n No seasons selected \n")
  if (max.age<min.age)      stop("\n Max age < min age   \n")

  cat("\n")
  cat("run.name         =",run.name,"\n")
  cat("n.persons        =",n.persons,"\n")
  cat("person.output    =",person.output,"\n")
  cat("source.output    =",source.output,"\n")
  cat("min.age          =",min.age,"\n")
  cat("max.age          =",max.age,"\n")
  cat("genders          =",genders,"\n")
  cat("season           =",seasons,"\n")
  cat("details          =",details,"\n")
  cat("age.match.pct    =",age.match.pct,"\n")
  cat("run.seed         =",run.seed,"\n")
  cat("set.size         =",set.size,"\n")
  cat("act.diary.file   =",act.diary.file,"\n")
  cat("chem.props.file  =",chem.props.file,"\n")
  cat("diet.diary.file  =",diet.diary.file,"\n")
  cat("exp.factor.file  =",exp.factor.file,"\n")
  cat("fugacity.file    =",fugacity.file,"\n")
  cat("media.file       =",media.file,"\n")
  cat("physiology.file  =",physiology.file,"\n")
  cat("population.file  =",population.file,"\n")
  cat("source.vars.file =",source.vars.file,"\n")
  cat("source.scen.file =",source.scen.file,"\n")
  cat("source.chem.file =",source.chem.file,"\n")
  cat("# chemicals      =",n.chem,"\n")

  s <- list(
    run.name         = run.name,
    n.persons        = n.persons,
    person.output    = person.output,
    source.output    = source.output,
    min.age          = min.age,
    max.age          = max.age,
    genders          = genders,
    seasons          = seasons,
    details          = details,
    age.match.pct    = age.match.pct,
    run.seed         = run.seed,
    set.size         = set.size,
    act.diary.file   = act.diary.file,
    chem.props.file  = chem.props.file,
    diet.diary.file  = diet.diary.file,
    exp.factor.file  = exp.factor.file,
    fugacity.file    = fugacity.file,
    media.file       = media.file,
    physiology.file  = physiology.file,
    population.file  = population.file,
    source.vars.file = source.vars.file,
    source.scen.file = source.scen.file,
    source.chem.file = source.chem.file,
    n.chem           = n.chem,
    chem.list        = chem.list)
  return(s)
}


#' read.act.diaries
#' 
#' Read.act.diaries reads human activity diaries from the .csv file indicated by "filename".  "Specs" contains the run specifications from the run.file, and is used to subset the activity diaries by age, gender, or season, if requested.  
#' 
#' @export
read.act.diaries = function(filename,specs) {
  # Read activity diaries file
  dt <- as.data.table(fread(paste0("inputs/",filename),colClasses = ("gender"="character")))
  if (nrow(dt)==0) stop("No data on diaries file \n")
  setnames(dt,names(dt),tolower(names(dt)))
  cols <- names(dt)[-which(names(dt) %in% c("chadid","gender","season","day.of.week"))]
  dt[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
  if(exists("veh.min",dt)) setnames(dt,"veh.min","in.veh.min")
  if(exists("veh.pai",dt)) setnames(dt,"veh.pai","in.veh.pai")
  d.min.age <- specs$min.age-round(specs$min.age*specs$age.match.pct/100)
  d.max.age <- specs$max.age+round(specs$max.age*specs$age.match.pct/100)
  if(d.max.age==0 & specs$age.match.pct>0) d.max.age <- 1
  mode(dt$age) <- "numeric"
  d <- dt[dt$age>=d.min.age & dt$age<=d.max.age & dt$season %in% specs$seasons
          & dt$gender %in% specs$genders]
    setnames(d,"age","d.age")
  setnames(d,"gender","d.gender")
  setnames(d,"weekend","d.weekend")
  setnames(d,"season","d.season")
  setnames(d,"bath","bath.mins")
  d[,diary.id:=1:nrow(d)]
  setkey(d,diary.id)
  cat("\n Reading Activity Diaries completed")
  return(d)
}


#' read.chem.props
#' 
#' Read.chem.props reads the chemical properties from the .csv file indicated by "filename". "Specs" contains the run specifications from the run.file, and is used to subset the chemicals by the list provided on the run.file.  If no such list was given, then all chemicals are kept.
#' 
#' @export
read.chem.props = function(filename,specs) {
  # Read chemical properties input file
  dt <- as.data.table(fread(paste0("inputs/",filename),na.strings=c("",".",".",".","NA")))
  setnames(dt,names(dt),tolower(names(dt)))
  if(exists("chemical",dt)) setnames(dt,"chemical","cas")
  mode(dt$cas) <- "character"
  mode(dt$kp) <- "numeric"
  dt$cas <- trimzero(str_replace_all(dt$cas,"-","_"))
  dt$cas <- ifelse(substr(dt$cas,1,1)=="_",paste0("0",dt$cas),dt$cas)
  if(specs$n.chem>0) dt <- dt[dt$cas %in% specs$chem.list]
  dt <- unique(dt,by="cas",fromLast=TRUE)
  if(nrow(dt)==0) cat("\n No chemicals left to model in source_chemicals \n")
  if(exists("name",dt))             setnames(dt,"name","chem.name")
  if(exists("vp.pa",dt))            setnames(dt,"vp.pa","vapor")
  if(exists("water.sol.mg.l",dt))   setnames(dt,"water.sol.mg.l","solub")
  if(exists("log.kow",dt))          dt$kow     <- 10^dt$log.kow
  if(exists("half.sediment.hr",dt)) dt$decay.s <- 24*log(2)/dt$half.sediment.hr
  if(exists("half.air.hr",dt))      dt$decay.a <- 24*log(2)/dt$half.air.hr
  cat("\n Reading Chemical Properties completed")
  return(dt)
}


#' read.diet.diaries
#' 
#' Read.diet.diaries reads the food diaries from the .csv file indicated by "filename". "Specs" contains the run specifications from the run.file, and is used to subset the activity diaries by age and/or gender.
#' 
#' @export
read.diet.diaries = function(filename,specs) {
  # Read diet diary consumption input file
  dt <- as.data.table(fread(paste0("inputs/",filename),header=TRUE))
  # setnames(dt,5:ncol(dt),toupper(names(dt)[5:ncol(dt)]))
  setnames(dt,names(dt),tolower(names(dt)))
  min.age     <- specs$min.age-round(specs$min.age*specs$age.match.pct/100)
  max.age     <- specs$max.age+round(specs$max.age*specs$age.match.pct/100)
  if(max.age==0 & specs$age.match.pct>0) max.age <- 1
  mode(dt$age) <- "numeric"
  dd           <- dt[dt$age>=min.age & dt$age<=max.age & dt$gender %in% specs$genders]
  setnames(dd,"age","f.age")
  setnames(dd,"gender","f.gender")
  dd[,diet.id:=1:nrow(dd)]
  setkey(dd,diet.id)
  cat("\n Reading Dietary Diaries completed")
  return(dd)
}


#' read.exp.factors
#' 
#' Read.exp.factors reads the exposure factors from the .csv file indicated by "filename".  
#' 
#' @export
read.exp.factors = function(filename) {
  # Read exposure factors input file
  df <- read.csv(paste0("inputs/",filename),as.is=TRUE)
  names(df) <- tolower(names(df))
  df <- df[!is.na(df$varname),]
  if (nrow(df)==0) stop ("No data on exposure factors file \n")
  df <- transform(df, form = tolower(gsub(" ", "", form)))
  mode(df$min.age) <- "numeric"
  mode(df$max.age) <- "numeric"
  df$min.age[is.na(df$min.age)] <- 0
  df$max.age[is.na(df$max.age)] <- 99
  mode(df$gender) <- "character"
  mode(df$season) <- "character"
  mode(df$media)  <- "character"
  df$gender <- toupper(as.character(df$gender))
  df$season <- toupper(as.character(df$season))
  df$media  <- tolower(as.character(df$media))
  dt <- as.data.table(df)
  setkey(dt, varname)
  dt$row <- 1:nrow(dt)
  cat("\n Reading Exposure Factors completed")
  return(dt)
}


#' read.fug.inputs
#' 
#' Read.fug.inputs reads the non-chemical dependent inputs for fugacity modeling in SHEDS from the .csv file indicated by "filename".   
#' 
#' @export
read.fug.inputs = function(filename) {
  # Read fugacity input file
  df <- read.csv(paste0("inputs/",filename),as.is=TRUE)
  names(df) <- tolower(names(df))
  if (nrow(df)==0) stop ("No data on fugacity input file \n")
  options(warn=-1)
  mode(df$par1) <- "numeric"
  mode(df$par2) <- "numeric"
  mode(df$par3) <- "numeric"
  mode(df$par4) <- "numeric"
  mode(df$lower.trun) <- "numeric"
  mode(df$upper.trun) <- "numeric"
  options(warn=0)
  df <- transform(df, form = tolower(gsub(" ", "", form)))
  dt <- as.data.table(df)
  dt$varname <- tolower(dt$varname)
  # cat("\n Reading Fugacity Variables completed")
  return(dt)
}


#' read.media.file
#' 
#' Read.media.file reads the names, properties, and associated microenvironments for each of the potential exposure media in SHEDS.
#' 
#' @export
read.media.file = function(filename) {
  # Read media input file
  dt <- as.data.table(read.csv(paste0("inputs/",filename),as.is=TRUE))
  setnames(dt,tolower(names(dt)))
  mode(dt$contact.p) <- "numeric"
  dt <- dt[dt$contact.p>0]
  if (nrow(dt)==0) stop("No data on exposure media file \n")
  dt$micro <- tolower(str_trim(dt$micro))
  dt$media <- tolower(str_trim(dt$media))
  dt$type  <- tolower(str_trim(dt$type))
  setkey(dt, micro)
  cat("\n Reading Media File completed")
  return(dt)
}



#' read.phys.file
#' 
#' Read.phys.file reads the physiology data from the .csv file indicated by "filename".
#' 
#' @export
read.phys.file = function(filename) {
  # Read physiology input file
  df <- read.csv(paste0("inputs/",filename),as.is=TRUE)
  names(df) <- tolower(names(df))
  mode(df$age) <- "numeric"
  dt <- as.data.table(df[!is.na(df$age),])
  if (nrow(dt)==0) stop("No data on physiology file \n")
  setkeyv(dt, c("gender","age"))
  cat("\n Reading Physiology File completed")
  return(dt)
}


#' read.pop.file
#' 
#' Read.pop.file reads the population data from the .csv file indicated by "filename".  
#' 
#' @export
read.pop.file = function(filename,specs) {
  # Read population input file
  df <- read.csv(paste0("inputs/",filename),as.is=TRUE)
  names(df) <- tolower(names(df))
  pop <- as.data.table(df[!is.na(df$age),])
  if (nrow(pop)==0) stop("No data on population input file \n")
  if (!"M" %in% specs$gender) pop$males   <- 0L
  if (!"F" %in% specs$gender) pop$females <- 0L
  mode(pop$age) <- "numeric"
  pop <- pop[pop$age>=specs$min.age & pop$age<=specs$max.age]
  if (sum(pop$males)+sum(pop$females)==0) stop("No population left\n")
  setkey(pop, age)
  cat("\n Reading Population File completed")
  return(pop)
}

#' read.source.scen.file
#' 
#' Read.source.scen.file reads the list of active exposure scenarios for each potential source of chemical from the .csv file indicated by "filename".  
#' 
#' @export
read.source.scen.file = function(filename) {
  # Read variables from srcScen input file
  fail <- FALSE
  df <- read.csv(paste0("inputs/",filename),as.is=TRUE)
  if (nrow(df)==0) stop ("No data on source_scenarios file \n")
  names(df)   <- tolower(str_trim(str_replace_all(names(df),","," ")))
  names(df)[substr(names(df),1,9)=="source.id"]  <- "src"
  names(df)[substr(names(df),1,8)=="food.res"]  <- "dietary"
  names(df)[substr(names(df),1,8)=="food.mig"]  <- "migration"
  names(df)[substr(names(df),1,21)=="product.direct.dermal"]  <- "dirderm"
  names(df)[substr(names(df),1,21)=="product.direct.ingest"]  <- "diringest"
  names(df)[substr(names(df),1,28)=="product.direct.inhalationaer"]  <- "dirinhaer"
  names(df)[substr(names(df),1,28)=="product.direct.inhalationvap"]  <- "dirinhvap"
  names(df)[substr(names(df),1,20)=="product.downthedrain"]  <- "downthedrain"
  names(df)[substr(names(df),1,16)=="product.indirect"]  <- "indir.fug"
  names(df)[substr(names(df),1,16)=="article.emission"]  <- "indir.y0"
  if(!exists("src",df))     {fail=TRUE; cat("\n No sources on srcScen file")}
  if(fail==TRUE) stop("\n Missing scenarios on source.scenarios file \n")
  df$src <- tolower(str_trim(str_replace_all(df$src,","," ")))
  df$source.type <- tolower(df$source.type)
  if(length(unique(df$src))!=nrow(df)) stop("\n Duplicate source.ID found")
  df <- check.src.scen.flags(df)
  df <- check.src.scen.types(df)
  dt <- as.data.table(df)
  dt$nscen <- dt$dietary + dt$dirderm + dt$diringest + dt$dirinhaer +
    dt$dirinhvap + dt$downthedrain + dt$indir.fug + dt$indir.y0 + dt$migration
  return(dt[dt$nscen>0])
}


#' check.src.scen.flags
#' 
#' This function checks whether the settings on the source.scen object (input as "df") are set to numeric values of 1 or 0.  If the column for a particular exposure scenario is missing from "df", it is created and assigned "0" for all sources.
#' 
#' @export
check.src.scen.flags = function(df) {
  fail <- FALSE
  if(mode(df$dietary)!="numeric") {
    if(mode(df$dietary)=="NULL") {df$dietary <- 0 }
    df$dietary <- tolower(str_trim(str_replace_all(df$dietary,","," ")))
    df$dietary[df$dietary=="." | df$dietary==""] <- 0
    mode(df$dietary) <- "numeric"
  }
  x <- length(unique(df$dietary))
  if((x!=1 && x!=2) || (x==2 && min(df$dietary)!=0 && max(df$dietary)!=1) ||
     (x==1 && (!any(df$dietary==0) && !any(df$dietary==1)))) {
    fail=TRUE; cat("\n Bad dietary flags") }

  if(mode(df$dirderm)!="numeric") {
    if(mode(df$dirderm)=="NULL") {df$dirderm <- 0 }
    df$dirderm <- tolower(str_trim(str_replace_all(df$dirderm,","," ")))
    df$dirderm[df$dirderm=="." | df$dirderm==""] <- 0
    mode(df$dirderm) <- "numeric"
  }
  x <- length(unique(df$dirderm))
  if((x!=1 && x!=2) || (x==2 && min(df$dirderm)!=0 && max(df$dirderm)!=1) ||
     (x==1 && (!any(df$dirderm==0) && !any(df$dirderm==1)))) {
    fail=TRUE; cat("\n Bad dirderm flags") }

  if(mode(df$diringest)!="numeric") {
    if(mode(df$diringest)=="NULL") {df$diringest <- 0 }
    df$diringest <- tolower(str_trim(str_replace_all(df$diringest,","," ")))
    df$diringest[df$diringest=="." | df$diringest==""] <- 0
    mode(df$diringest) <- "numeric"
  }
  x <- length(unique(df$diringest))
  if((x!=1 && x!=2) || (x==2 && min(df$diringest)!=0 && max(df$diringest)!=1) ||
     (x==1 && (!any(df$diringest==0) && !any(df$diringest==1)))) {
    fail=TRUE; cat("\n Bad diringest flags") }

  if(mode(df$dirinhaer)!="numeric") {
    if(mode(df$dirinhaer)=="NULL") {df$dirinhaer <- 0 }
    df$dirinhaer <- tolower(str_trim(str_replace_all(df$dirinhaer,","," ")))
    df$dirinhaer[df$dirinhaer=="." | df$dirinhaer==""] <- 0
    mode(df$dirinhaer) <- "numeric"
  }
  x <- length(unique(df$dirinhaer))
  if((x!=1 && x!=2) || (x==2 && min(df$dirinhaer)!=0 && max(df$dirinhaer)!=1) ||
     (x==1 && (!any(df$dirinhaer==0) && !any(df$dirinhaer==1)))) {
    fail=TRUE; cat("\n Bad dirinhaer flags") }

  if(mode(df$dirinhvap)!="numeric") {
    if(mode(df$dirinhvap)=="NULL") {df$dirinhvap <- 0 }
    df$dirinhvap <- tolower(str_trim(str_replace_all(df$dirinhvap,","," ")))
    df$dirinhvap[df$dirinhvap=="." | df$dirinhvap==""] <- 0
    mode(df$dirinhvap) <- "numeric"
  }
  x <- length(unique(df$dirinhvap))
  if((x!=1 && x!=2) || (x==2 && min(df$dirinhvap)!=0 && max(df$dirinhvap)!=1) ||
     (x==1 && (!any(df$dirinhvap==0) && !any(df$dirinhvap==1)))) {
    fail=TRUE; cat("\n Bad dirinhvap flags") }

  if(mode(df$downthedrain)!="numeric") {
    if(mode(df$downthedrain)=="NULL") {df$downthedrain <- 0 }
    df$downthedrain <- tolower(str_trim(str_replace_all(df$downthedrain,","," ")))
    df$downthedrain[df$downthedrain=="." | df$downthedrain==""] <- 0
    mode(df$downthedrain) <- "numeric"
  }
  x <- length(unique(df$downthedrain))
  if((x!=1 & x!=2) || (x==2 & min(df$downthedrain)!=0 & max(df$downthedrain)!=1)
     || (x==1 && (!any(df$downthedrain==0) && !any(df$downthedrain==1)))) {
    fail=TRUE; cat("\n Bad downthedrain flags") }

  if(mode(df$indir.fug)!="numeric") {
    if(mode(df$indir.fug)=="NULL") {df$indir.fug <- 0 }
    df$indir.fug <- tolower(str_trim(str_replace_all(df$indir.fug,","," ")))
    df$indir.fug[df$indir.fug=="." | df$indir.fug==""] <- 0
    mode(df$indir.fug) <- "numeric"
  }
  x <- length(unique(df$indir.fug))
  if((x!=1 & x!=2) || (x==2 & min(df$indir.fug)!=0 && max(df$indir.fug)!=1) ||
     (x==1 && (!any(df$indir.fug==0) && !any(df$indir.fug==1)))) {
    fail=TRUE; cat("\n Bad indir.fug flags") }

  if(mode(df$indir.y0)!="numeric") {
    if(mode(df$indir.y0)=="NULL") {df$indir.y0 <- 0 }
    df$indir.y0 <- tolower(str_trim(str_replace_all(df$indir.y0,","," ")))
    df$indir.y0[df$indir.y0=="." | df$indir.y0==""] <- 0
    mode(df$indir.y0) <- "numeric"
  }
  x <- length(unique(df$indir.y0))
  if((x!=1 & x!=2) | (x==2 && min(df$indir.y0)!=0 && max(df$indir.y0)!=1) |
     (x==1 & (!any(df$indir.y0==0) & !any(df$indir.y0==1)))) {
    fail=TRUE; cat("\n Bad indir.y0 flags") }

  if(fail==TRUE) stop("\n Bad option flags on source.scenarios file \n")
  return(df)
}


#' Check.src.scen.types
#' 
#' This function checks whether the settings on the source.scen object (input as "df") are consistent with the source.type for each source.
#' 
#' @export
check.src.scen.types = function(df) {
  if(any(df$dietary[!df$source.type=="food"]==1 )) {
    df$dietary[!df$source.type=="food"] <- 0
    cat("\n Non-food types flagged for food.residue scenario \n")
  }
  if(any(df$migration[!df$source.type=="food"]==1 )) {
    df$migration[!df$source.type=="food"] <- 0
    cat("\n Non-food types flagged for food.migration scenario \n")
  }
  if(any(df$dirderm[!df$source.type=="product"]==1)) {
    df$dirderm[!df$source.type=="product"] <- 0
    cat("\n Non-product types flagged for product.direct.dermal scenario \n")
  }
  if(any(df$diringest[!df$source.type=="product"]==1)) {
    df$diringest[!df$source.type=="product"] <- 0
    cat("\n Non-product types flagged for product.direct.ingestion scenario \n")
  }
  if(any(df$dirinhaer[!df$source.type=="product"]==1)) {
    df$dirinhaer[!df$source.type=="product"] <- 0
    cat("\n Non-product types flagged for product.direct.inhalationaerosol scenario \n")
  }
  if(any(df$dirinhvap[!df$source.type=="product"]==1)) {
    df$dirinhvap[!df$source.type=="product"] <- 0
    cat("\n Non-product types flagged for product.direct.inhalationvapor scenario \n")
  }
  if(any(df$downthedrain[!df$source.type=="product"]==1)) {
    df$downthedrain[!df$source.type=="product"] <- 0
    cat("\n Non-product types flagged for product.downthedrain scenario \n")
  }
  if(any(df$indir.fug[!df$source.type=="product"]==1)) {
    df$indir.fug[!df$source.type=="product"] <- 0
    cat("\n Non-product types flagged for product.indirect scenario \n")
  }
  if(any(df$indir.y0[!df$source.type=="article"]==1)) {
    df$indir.y0[!df$source.type=="article"] <- 0
    cat("\n Non-article types flagged for article.emission scenario \n")
  }
  if(any(df$indir.fug[df$indoor==0]==1)) {
    df$indir.fug[df$indoor==0] <- 0
    cat("\n Non-indoor source flagged for product.indirect scenario \n")
  }
  if(any(df$indir.y0[df$indoor==0]==1)) {
    df$indir.y0[df$indoor==0] <- 0
    cat("\n Non-indoor source flagged for article.emission scenario \n")
  }
  return(df)
}


#' read.source.chem.file
#' 
#' Read.source.chem.file reads the distributions that are specific to combinations of source and chemical from the indicated .csv file. 
#' 
#' @export
read.source.chem.file = function(filename,scenSrc,specs) {
  # Read variables from srcChem input file
  fail <- FALSE
  df <- read.csv(paste0("inputs/",filename),as.is=TRUE)
  test2<<-df
  if (nrow(df)==0) stop ("No data on source.chemicals file \n")
  names(df)   <- tolower(str_trim(str_replace_all(names(df),","," ")))
  names(df)[substr(names(df),1,9)=="source.id"]  <- "src"
  names(df)[substr(names(df),1,9)=="source.de"]  <- "description"
  names(df)[substr(names(df),1,3)=="cas"]  <- "cas"
  names(df)[substr(names(df),1,3)=="che"]  <- "cas"
  names(df)[substr(names(df),1,3)=="var"]  <- "varname"
  names(df)[substr(names(df),1,3)=="des"]  <- "description"
  names(df)[substr(names(df),1,3)=="uni"]  <- "units"
  names(df)[substr(names(df),1,3)=="gen"]  <- "gender"
  names(df)[substr(names(df),1,3)=="min"]  <- "min.age"
  names(df)[substr(names(df),1,3)=="max"]  <- "max.age"
  names(df)[substr(names(df),1,3)=="low"]  <- "lower.trun"
  names(df)[substr(names(df),1,3)=="upp"]  <- "upper.trun"
  names(df)[substr(names(df),1,3)=="for"]  <- "form"
  names(df)[substr(names(df),1,3)=="mea"]  <- "mean"
  names(df)[substr(names(df),1,3)=="val"]  <- "values"
  names(df)[substr(names(df),1,2)=="cv"]   <- "cv"
  if(length(names(df))>length(unique(names(df)))) {
    fail <- TRUE
    cat("\n Duplicate variable names on sources file")
  }
  if(!exists("src",df))     {fail<-TRUE; cat("\n No source.ID on src.chem file")}
  if(!exists("cas",df))     {fail<=TRUE; cat("\n No cas numbers on src.chem file")}
  if(!exists("varname",df)) {fail<-TRUE; cat("\n No varname on src.chem file")}
  if(!exists("form",df))    {fail<-TRUE; cat("\n No form on src.chem file")}
  df$form <- tolower(str_trim(str_replace_all(df$form,","," ")))
  if(any(df$form=="empirical") && !exists("values",df)) {
    fail <- TRUE
    cat("\n Values not given with empirical distribution")
  }
  if(any(df$form!="empirical") && (!exists("mean",df) | !exists("cv",df))) {
    fail <- TRUE
    cat("\n Mean and CV not given with non-empirical distribution")
  }
  if(!exists("description",df)) df$description <- ""
  if(!exists("units",df))       df$units       <- ""
  if(!exists("gender",df))      df$gender      <- ""
  if(!exists("min.age",df))     df$min.age     <- 0
  if(!exists("max.age",df))     df$max.age     <- 99
  if(!exists("lower.trun",df))  df$lower.trun  <- 0
  if(!exists("upper.trun",df))  df$upper.trun  <- 1E6
  if(!exists("mean",df))        df$mean        <- 0
  if(!exists("cv",df))          df$cv          <- 0
  if(!exists("values",df))      df$values      <- ""
  df$src      <- tolower(str_trim(str_replace_all(df$src,","," ")))
  df$cas      <- tolower(str_trim(str_replace_all(df$cas,","," ")))
  df$varname  <- tolower(str_trim(str_replace_all(df$varname,","," ")))
  df$descrip  <- tolower(str_trim(str_replace_all(df$description,","," ")))
  df$units    <- tolower(str_trim(str_replace_all(df$units,","," ")))
  df$gender   <- toupper(str_trim(str_replace_all(df$gender,","," ")))
  mode(df$min.age) <- "numeric"
  mode(df$max.age) <- "numeric"
  df$min.age[is.na(df$min.age)] <- 0
  df$max.age[is.na(df$max.age)] <- 99
  if(any(df$min.age<0)) {fail<-TRUE; cat("\n Negative age on sources file") }
  df$max.age[df$max.age>99]     <- 99
  df$gender[df$gender=="W"] <- "F"
  df$gender[df$gender=="B"] <- ""
  df$gender[df$gender==" "] <- ""
  df <- df[!is.na(df$cas),]
  df$cas <- trimzero(str_replace_all(df$cas,"-","_"))
  if(mode(df$mean)!="numeric") {
    df$mean <- tolower(str_trim(str_replace_all(df$mean,","," ")))
    df$mean[df$mean=="." | df$mean==""] <- -1
    mode(df$mean) <- "numeric"
    df$mean[df$mean<0] <- NA
  }
  if(mode(df$cv)!="numeric") {
    df$cv <- tolower(str_trim(str_replace_all(df$cv,","," ")))
    df$cv[df$cv=="." | df$cv==""] <- -1
    mode(df$cv) <- "numeric"
    df$cv[df$cv<0] <- NA
  }
  if(any(df$mean<0,na.rm=TRUE)) {fail<-TRUE; cat("\n Negative mean on srcChem file")}
  if(any(df$cv<0,na.rm=TRUE))   {fail<-TRUE; cat("\n Negative cv on srcChem file") }
  dt  <- as.data.table(df)
  if (specs$n.chem>0) { dt <- dt[dt$cas %in% specs$chem.list]
  } else {dt <- dt[dt$cas %in% chem.props$cas]}
  chem.list  <- sort(unique(dt$cas))
  dt <- dt[dt$src %in% scenSrc]
  dt[substr(dt$varname,1,2)=="f."]$upper.trun <- 1
  dt$form[dt$form=='lognormal' & is.na(dt$cv)] <- 'point'
  dt$form[dt$form=='normal' & dt$cv>0.5] <- 'lognormal'
  dt$values[dt$form=='empirical' & is.na(dt$values)] <- "1"
  if (fail==TRUE) stop()
  cat("\n Reading Source.chemicals file completed")
  return(dt)
}


#' update.specs
#' 
#' "Specs" is the list of run settings read from the run.file.  "Dt" is a data table of source-chemical combinations for which distributions have been specified.  "Specs" contains a list of chemicals to be processed in the SHEDS run, and if any of these chemicals are missing from the "dt" table, then update.specs removes them from the list.  Otherwise, "specs" is not altered.
#' 
#' @export
update.specs = function(specs,dt) {
  specs$chem.list <- sort(unique(dt$cas))
  specs$n.chem    <- length(specs$chem.list)
  if (specs$n.chem==0) {cat ("\n No chemicals remaining to model \n"); stop() }
  return(specs)
}


#' read.source.vars.file
#' 
#' Read.source.vars.file reads the distributions that are specific to combinations of source and chemical from the indicated .csv file.  
#' 
#' @export
read.source.vars.file = function(filename,src.scen) {
  # Read variables from source.variables input file
  fail <- FALSE
  df <- read.csv(paste0("inputs/",filename),as.is=TRUE)
  if (nrow(df)==0) stop ("No data on source.variables file \n")
  names(df)   <- tolower(str_trim(str_replace_all(names(df),","," ")))
  names(df)[substr(names(df),1,9)=="source.id"]  <- "src"
  names(df)[substr(names(df),1,9)=="source.de"]  <- "src.description"
  names(df)[substr(names(df),1,8)=="variable"]   <- "varname"
  names(df)[substr(names(df),1,8)=="var.desc"]   <- "var.description"
  names(df)[substr(names(df),1,3)=="uni"]  <- "units"
  names(df)[substr(names(df),1,3)=="gen"]  <- "gender"
  names(df)[substr(names(df),1,3)=="min"]  <- "min.age"
  names(df)[substr(names(df),1,3)=="max"]  <- "max.age"
  names(df)[substr(names(df),1,3)=="low"]  <- "lower.trun"
  names(df)[substr(names(df),1,3)=="upp"]  <- "upper.trun"
  names(df)[substr(names(df),1,3)=="res"]  <- "resamp"
  names(df)[substr(names(df),1,3)=="pro"]  <- "probs"
  names(df)[substr(names(df),1,3)=="val"]  <- "values"
  names(df)[substr(names(df),1,3)=="for"]  <- "form"
  names(df)[substr(names(df),1,3)=="mea"]  <- "mean"
  names(df)[substr(names(df),1,2)=="cv"]   <- "cv"
  if(length(names(df))>length(unique(names(df)))) {
    fail <- TRUE
    cat("\n Duplicate variable names on source.variables file")
  }
  if(!exists("src",df))     {fail<-TRUE; cat("\n No source on src.vars file")}
  if(!exists("varname",df)) {fail<-TRUE; cat("\n No varname on src.vars file")}
  if(!exists("form",df))    {fail<-TRUE; cat("\n No form on src.vars file")}
  df$form <- tolower(str_trim(str_replace_all(df$form,","," ")))
  if(any(df$form=="empirical") && !exists("values",df)) {
    fail <- TRUE
    cat("\n Values not given with empirical distribution")
  }
  if(any(df$form!="empirical") && (!exists("mean",df) | !exists("cv",df))) {
    fail <- TRUE
    cat("\n Mean and CV not given with non-empirical distribution")
  }
  if(!exists("description",df)) df$description <- ""
  if(!exists("units",df))       df$units       <- ""
  if(!exists("gender",df))      df$gender      <- ""
  if(!exists("min.age",df))     df$min.age     <- 0
  if(!exists("max.age",df))     df$max.age     <- 99
  if(!exists("lower.trun",df))  df$lower.trun  <- 0
  if(!exists("upper.trun",df))  df$upper.trun  <- 1E6
  if(!exists("mean",df))        df$mean        <- 0
  if(!exists("cv",df))          df$cv          <- 0
  if(!exists("values",df))      df$values      <- ""
  if(!exists("probs",df))       df$probs       <- 0
  if(!exists("resamp",df))      df$resamp      <- "y"
  df$src      <- tolower(str_trim(str_replace_all(df$src,","," ")))
  df$varname  <- tolower(str_trim(str_replace_all(df$varname,","," ")))
  df$descrip  <- tolower(str_trim(str_replace_all(df$description,","," ")))
  df$units    <- tolower(str_trim(str_replace_all(df$units,","," ")))
  df$gender   <- toupper(str_trim(str_replace_all(df$gender,","," ")))
  mode(df$min.age) <- "numeric"
  mode(df$max.age) <- "numeric"
  df$min.age[is.na(df$min.age)] <- 0
  df$max.age[is.na(df$max.age)] <- 99
  if(any(df$min.age<0)) {fail<-TRUE; cat("\n Negative age on srcVars file") }
  df$max.age[df$max.age>99]     <- 99
  df$gender[df$gender=="W"] <- "F"
  df$gender[df$gender=="B"] <- ""
  df$gender[df$gender==" "] <- ""
  if(mode(df$mean)!="numeric") {
    df$mean <- tolower(str_trim(str_replace_all(df$mean,","," ")))
    df$mean[df$mean=="." | df$mean==""] <- -1
    mode(df$mean) <- "numeric"
    df$mean[df$mean<0] <- NA
  }
  if(mode(df$cv)!="numeric") {
    df$cv <- tolower(str_trim(str_replace_all(df$cv,","," ")))
    df$cv[df$cv=="." | df$cv==""] <- -1
    mode(df$cv) <- "numeric"
    df$cv[df$cv<0] <- NA
  }
  if(any(df$mean<0,na.rm=TRUE)) {fail<-TRUE; cat("\n Negative mean on srcVars file")}
  if(any(df$cv<0,na.rm=TRUE))   {fail<-TRUE; cat("\n Negative cv on srcVars file") }

  dt <- as.data.table(df)
  dt <- dt[dt$src %in% src.scen$src]
  dt[substr(dt$varname,1,2)=="f."]$upper.trun <- 1
  dt$form[dt$form=='lognormal' & is.na(dt$cv)] <- 'point'
  dt$form[dt$form=='normal' & dt$cv>0.5] <- 'lognormal'
  for (i in 1:nrow(src.scen)) {
    sc <- src.scen[i]
    v1 <- dt[dt$src==sc$src]
    if (any(substr(v1$varname,1,8)=="use.prev"))  u.prev<-TRUE else u.prev<-FALSE
    if (any(substr(v1$varname,1,8)=="use.freq"))  u.freq<-TRUE else u.freq<-FALSE
    if (any(substr(v1$varname,1,9)=="home.prev")) h.prev<-TRUE else h.prev<-FALSE
    if (any(substr(v1$varname,1,4)=="mass"))      mass  <-TRUE else mass  <-FALSE
    if (any(substr(v1$varname,1,3)=="dur"))       dur   <-TRUE else dur   <-FALSE
    if (any(substr(v1$varname,1,3)=="vol"))       vol   <-TRUE else vol   <-FALSE
    if (any(substr(v1$varname,1,4)=="area"))      area  <-TRUE else area  <-FALSE
    if (any(substr(v1$varname,1,6)=="f.area"))    f.area<-TRUE else f.area<-FALSE
    if (any(substr(v1$varname,1,2)=="y0"))        y0    <-TRUE else y0    <-FALSE
    if (any(substr(v1$varname,1,5)=="f.con"))     f.con <-TRUE else f.con <-FALSE
    if (any(substr(v1$varname,1,5)=="f.res"))     f.res <-TRUE else f.res <-FALSE
    if (any(substr(v1$varname,1,5)=="f.ing"))     f.ing <-TRUE else f.ing <-FALSE
    if (any(substr(v1$varname,1,5)=="f.aer"))     f.aer <-TRUE else f.aer <-FALSE
    if (any(substr(v1$varname,1,5)=="f.dra"))     f.dra <-TRUE else f.dra <-FALSE
    if (sc$dirderm==1) {
      if (u.prev==FALSE) cat ("\n Use.prev missing for source ", sc$src)
      if (u.freq==FALSE) cat ("\n Use.freq missing for source ", sc$src)
      if (mass  ==FALSE) cat ("\n Mass missing for source ", sc$src)
      if (f.con ==FALSE) cat ("\n F.contact missing for source ", sc$src)
      if (f.res ==FALSE) cat ("\n F.residual missing for source ", sc$src)
    }
    if (sc$diringest==1) {
      if (u.prev==FALSE) cat ("\n Use.prev missing for source ", sc$src)
      if (u.freq==FALSE) cat ("\n Use.freq missing for source ", sc$src)
      if (mass  ==FALSE) cat ("\n Mass missing for source ", sc$src)
      if (f.ing ==FALSE) cat ("\n F.ingested missing for source ", sc$src)
    }
    if (sc$dirinhaer==1) {
      if (u.prev==FALSE) cat ("\n Use.prev missing for source ", sc$src)
      if (u.freq==FALSE) cat ("\n Use.freq missing for source ", sc$src)
      if (mass  ==FALSE) cat ("\n Mass missing for source ", sc$src)
      if (dur   ==FALSE) cat ("\n Duration missing for source ", sc$src)
      if (vol   ==FALSE) cat ("\n Volume missing for source ", sc$src)
      if (f.aer ==FALSE) cat ("\n F.aerosol missing for source ", sc$src)
    }
    if (sc$dirinhvap==1) {
      if (u.prev==FALSE) cat ("\n Use.prev missing for source ", sc$src)
      if (u.freq==FALSE) cat ("\n Use.freq missing for source ", sc$src)
      if (mass  ==FALSE) cat ("\n Mass missing for source ", sc$src)
      if (dur   ==FALSE) cat ("\n Duration missing for source ", sc$src)
    }
    if (sc$downthedrain==1) {
      if (u.prev==FALSE) cat ("\n Use.prev missing for source ", sc$src)
      if (u.freq==FALSE) cat ("\n Use.freq missing for source ", sc$src)
      if (mass  ==FALSE) cat ("\n Mass missing for source ", sc$src)
      if (f.dra ==FALSE) cat ("\n F.drain missing for source ", sc$src)
    }
    if (sc$indir.fug==1) {
      if (h.prev==FALSE & u.prev==FALSE)
        cat ("\n Home.prev and use.prev missing for source", sc$src)
      if (u.freq==FALSE) cat ("\n Use.freq missing for source ", sc$src)
      if (mass  ==FALSE) cat ("\n Mass missing for source ", sc$src)
    }
    if (sc$indir.y0==1) {
      if (h.prev==FALSE & u.prev==FALSE)
        cat ("\n Home.prev and use.prev missing for source ", sc$src)
      if (f.area==FALSE) cat ("\n Area missing for source ", sc$src)
    }
  }
  if (fail==TRUE) stop("\n Bad data on source.variables file \n")
  cat("\n Reading Source.variables file completed")
  tester<<-dt
  return(dt)
}


#' check.foods
#'  
#' This function creates a list of unique food types.
#' 
#' @export
check.foods = function(s) {
  t <- unique(s,fromLast=TRUE,by="varname")
  return(t)
}


# Notes for SHEDS-HT ReadData module
# Last modified by WGG on June 4, 2015
#
# read.run.file()   Apart from running steup(), this is the first step in any
#                 SHEDS run. All of the run-specific information is on the
#                 "run file", which is a text file that must be in the SHEDS
#                 /inputs folder.  If no file name is specified in the call,
#                 SHEDS uses the generic "runFile.csv".  All input files must
#                 have the .csv file extension, but this may be omitted from
#                 the function call, so run("run.file") is equivalent to
#                 run("runFile.csv"). Users may want to construct many run
#                 files for specific purposes; for example, to examine subsets
#                 of chemicals, or scenarios, or of the population, like
#                 particular age groups. Natually, each such file must have
#                 a unique name, which is suppled as an argument to run().
#
#                 The run file contains settings for 12 model parameters and
#                 9 input file names.  Defaults are provided, but it is good
#                 practice to specify all values in the run file, even if
#                 the default is used. The first line of the run file is a
#                 header, telling SHEDS the names of the two data columns. An
#                 equal sign "=" is used as the separator between the variable
#                 name and its setting.
#
#                 For inputs such as gender and season, the input is read as a
#                 string which is then searched for key letters. If "f" or "F"
#                 is found in the gender string, then females are run. If "m"
#                 or "M" is found, then males are run.  For seasons, w=winter,
#                 p=spring, s=summer, and f=fall or autumn.
#
#                 The output from read.run.file() is the global object "specs"
#                 which contains 21 components.  This is the only global item
#                 in a normal SHEDS run. The user may choose to create others
#                 to aid in debugging. All SHEDS functions use explicitly
#                 passed arguments for all inputs and outputs, except for
#                 references to "specs".
#
# read.act.diaries()  This function reads the activity diary database. The first
#                   line of the file contains the column headers (variable
#                   names), which should be chadid, age, gender, day.of.week,
#                   weekend, season, month, bath, InAwkMin, InSlpMin, InWrkMin,
#                   InVehMin, InOthMin, OutHmMin, OutWrkMin, OutOthMin,
#                   InAwkPai, InSlpPai, InWrkPai, InVehPai, InOthPai, OutHmPai,
#                   OutWrkPai, and OutOthPai. The 7 variables ending in "Min"
#                   are the times in minutes spent in each of the 7 micros,
#                   and together must sum to 1440 minutes.  The 7 ending in
#                   "Pai" are the physical activity indices for those activities.
#                   The chadid may be any character string. All data for each
#                   person appear on one line. Technically, the chadid does not
#                   have to be unique, and duplicates will appear as separate
#                   diaries. All acceptable diaries are assigned unique numbers.
#                   The age should be numeric, integer, from 0-99 inclusive.
#                   The gender should be M or F. The day.of.week is not
#                   currently used, but is meant to be a 3-character label
#                   such as SUN, MON, etc. The weekend is either 1=SAT or SUN,
#                   or 0=other day.  The season is W, P, S, or F. The month
#                   is numbered 1-12, but is not directly used by SHEDS.  The
#                   bath variable indicates the number of minutes reported as
#                   bath/shower activity.  For many diaries this is not
#                   separately reported and therefore is zero.  For a fraction
#                   of these zeroes, a bath/shower is randomly generated later
#                   on.
#
#                   Some diaries may be deleted from the database. The ones kept
#                   must match one of the modeled genders and seasons, and must
#                   be within the age.match.pct of the modeled ages.  For example,
#                   if specs$min.age=20, specs$max.age=40, and specs$age.match.pct=25,
#                   then diaries between ages 15 and 50 would be kept (that is,
#                   20*0.75 to 40*1.25).  These diary ages could be selected
#                   to represent the target population.  Diaries with no
#                   selection probability are dropped.
#
# read.chem.props()   This function reads the chemical properties file. Each
#                   chemical is identified by a CAS number.  If the CAS number
#                   contains hyphens, they are changed to underscores to avoid
#                   problems when opening such files in Excel.  If the same
#                   CAS number appears on multiple lines, only the last one
#                   is kept.
#
#                   CAS numbers may include leading zeroes, but these are
#                   sometimes removed, which leads to difficulties in matching
#                   numbers across files.  This function trims leading zeroes,
#                   except a single zero is retained if all the digits before
#                   the first underscore happen to be zero.
#
#                   Some variables, including the vapor pressure and solubility,
#                   are renamed. THe log.Kow value is converted to Kow.  The
#                   air and sediment decay rates are converted from half-lives
#                   in hours to daily decay rates.
#
# read.diet.diaries() Eacgh diet diary gives the consumption in grams per day for
#                   each food group, for the given person.  The person is
#                   identified by ID, gender, age, and body weight in kg. There
#                   are no variables for season or weekend.  The food groups
#                   are not predefined: SHEDS uses all the column headers for
#                   the fifth column and beyond as the names of the food groups.
#                   All are converted to uppercase, and cannot be longer than
#                   7 characters.  These are later paired with the food groups
#                   found on the scenarios input file.  As with the activity
#                   diaries, the food diaries that are possible matches to the
#                   age and gender of simulated people are kept, while other
#                   diaries are discarded.
#
# read.exp.factors()  The ExpFactors file contains definitions of distributions
#                   for various exposure factors such as transfer coefficients,
#                   hand-to-mouth variables, and dermal removal terms. Each one
#                   is defined by a name, eight parameters, and five demograpahic
#                   factors. The varname is the variable name as used in the
#                   SHEDScode, so that cannot be altered without changing the
#                   code as well. The eight parameters are form (or type of
#                   distribution), 4 numeric parameters, optional minimum and
#                   maximum bounds, and the truncation option. If this last
#                   option is "yes" or missing, then the distribution is
#                   effectively resampled until a random value between the
#                   bounds is selected.  Otherwise, values outside the bounds
#                   are moved to the appropriate bounding value.
#
#                   The five demographic variables are min.age, max.age, gender,
#                   season, and media. If any or all of these are missing, the
#                   distribution is assumed to apply to all suitable cases. Due
#                   to limited data, it is uncommon to use any restrictions other
#                   than age ranges.  The hand-to-mouth variables are almost
#                   always age-dependent.
#
# read.fug.inputs()   This function reads the file containing inputs for fugacity
#                   modeling variables.  These are similar to the ExpFactors, but
#                   the fugacity variables do not have demographic factors. The
#                   file therefore consists of lines with a varname, eight
#                   parameters, and an optional description at the end.
#
# read.media.file()   Each row on the media file has six columns: micro, micro.label,
#                   media, type, contact.p, and media.label.  The micro is the
#                   activity diary location where the medium is found. The person
#                   must spend time in that micro or else no exposure there can
#                   occur.  THe focus of SHEDS has initially been on the in-home
#                   micro, but others are included and may be used for exposure.
#                   The micro.label is not used by SHEDs. The type indicates the
#                   relevant exposure pathway: the current types are surface, air,
#                   and (not yet coded) pet.  Each combination of micro and type
#                   has a specific media name.
#
#                   Contact.p gives the fraction of time in that micro that could
#                   reasonably be expected to be potential contact time. For air
#                   and for surfaces found everywhere, contact.p should be one.
#                   For pets, and for specific surfaces (such as wooden decks
#                   outdoors), contact.p should be lower, reflecting the fact that
#                   much of the time in the relevant micro is spent away from that
#                   particular medium. The diary times are multiplied by contact.p
#                   before any further calculations are performed.  Setting contact.p
#                   to zero is the simplest way to limit SHEDS to a particular subset
#                   of all the possible media.
#
# read.phys.file()    The physiology input file contains regression parameters for
#                   body weight, height, and basal metabolic rate (BMRI), for each
#                   age and gender.  For body weight, the parameters apply to the
#                   natural logarithm of the weight (in kg), and the mean and
#                   standard deviation are given. For height, two approaches are
#                   used, depending on the person's age. Below age 20, height is
#                   is normally distributed around an age-gender specific mean.
#                   For people age 20 or older, the height is given by a regression
#                   on body-weight, also specific to age and gender. All input lines
#                   on this file allow space for the height mean, standard deviation,
#                   slope, intercept, and residual (five values in total), even
#                   though no line contains all five values. After the height
#                   parameters come three BMR parameters (slope, intercept, and
#                   residual).  The BMR slope is multiplied by the body weight.
#
# read.pop.file()   This function reads the population file.  This is a very simple
#                 file with age, #males, and #females on each line. These numbers
#                 are total counts from the 2000 U.S. Census.  The purpose of this
#                 is to create appropriate numbers of simulated persons of each age.
#                 For example, in a run from age 0-99, the number of 99-year old
#                 persons should be quite small compared to younger ages. Less
#                 extreme effects exist at other ages.
#
#                 The user may limit the run by age, gender, or both. In that case,
#                 the "counts" of the unallowed combinations are set to zero. Both
#                 age and gender are selected at random in SHEdS, but if the weight
#                 is zero, then that option cannot occur.
#
# readSrcScenFile() This is one of three input files defining chemical sources. Each row
#                   of the SrcScen file contains a "source", which may be a product,
#                   article, or food, and a series of switches. There is one switch for
#                   each exposure scenario in SHEDS.  The switch is set to 1 to indicate
#                   that scenarios exists for that source, or set to 0 otherwise. The
#                   input file may contain all switches set to 0 for some source, in
#                   which case SHEDS deletes that source from its list. Any combination
#                   of switches may be set to 1, even all of them for the same source.
#
#                   Editing the SrcScen file is the simplest way to control which sources
#                   and which exposure scenarios are active in a run. However, this file
#                   does not provide a way to limit which chemicals are modeled. The main
#                   way to control the chemicals is on the Run File.
#
# read.src.chem.file() This is the second of the input files defining chemical sources. This
#                   file contains distributions for variables that depend on both the source
#                   and on the chemical. For non-foods, these are the variables f.chemical
#                   and chem.prev. For food residues, these variables are detects, non-detects,
#                   and residue.
#
#                   Any combinations of source and chemical that are not listed on this file
#                   are assumed to be zero (that is, the given chemical is not found in that
#                   source).Two important lines in this function are
#
#                          if (specs$n.chem>0) dt <- dt[dt$cas %in% specs$chem.list]
#                   and
#                          dt <- dt[dt$src %in% scenSrc]
#
#                   The first of these limits the data to chemicals found in specs$chem.list,
#                   unless that list was empty (no chemicals listed on the Run File). The
#                   second line limits the data to sources on the SrcScen file.
#                   This allows the user to define a very large (i.e., universal) SrcChem file
#                   that can be used in all runs, while each SHEDS run extracts the relevant
#                   information from it.
#
# read.src.vars.file() This is like the read.src.chem.file function, except that it contains
#                   distributions for variables that depend only on the source and not on
#                   the chemical.  These variables include use.prev and home.prev which
#                   indicate market penetration, use.freq, mass, area, volume, duration,
#                   f.contact, f.residual, f.aerosol, and f.ingested, which describe product
#                   use.
#
# check.foods()  This function lists the food groups found under the current
#               Dietary scenario.


