#' combine_output
#' 
#' combine_output is a post-processing tool that extracts selected statistics (those specified in the 'metrics' argument)
#' from the 'allstats' file for each chemical in a single model run, and combines them onto out.file.   
#'  
#' @param run.name default = name of last run (in current R session)
#' 
#' @param out.file default="SHEDSOutFCS.csv"
#' 
#' @param metrics the summary statsistics to be pulled. default = c("5\%", "50\%", "75\%", "95\%", "99\%","mean","sd")
#' 
#' @return finaldata A R object that has "exp.dermal", "exp.ingest", m"exp.inhal", "dose.inhal", "dose.intake", "abs.dermal.ug",
#'         "abs.ingest.ug", "abs.inhal.ug", "abs.tot.ug", "abs.tot.mgkg", "ddd.mass" for the chemcial of question. This is tablated for 
#'         the "5\%", "50\%", "75\%", "95\%", "99\%", quantile of exposure as well as the mean and sd 
#'         \cr
#'         "finaldata,paste0("output/",run.name,"/",out.file.csv") this is the CSV of finaldata
#'         
#' @details In a multichemical run, SHEDS creates separate output files for each chemical.  If the analyst wants to compare
#'          exposure or dose metrics across chemicals, use combine_output to put the relevant information for all chemicals together
#'          on one file.
#' 
#' @author  Kristin Isaacs, Graham Glen
#'   
#' @export

# Functions for postprocessing SHEDS-HT output files
# All functions written by KKI
# Last changes by WGG Aug 28, 2016

combine_output = function(run.name=specs$run.name,out.file="SHEDSOutFCS.csv",metrics=c("5%","50%","75%","95%","99%","mean","sd") ) {
  
  allfiles   <-list.files(path=paste0("output/",run.name), all.files=FALSE, include.dirs=FALSE, full.names=TRUE)
  allcasnums <- list()
  keep  <- rep(FALSE,length(allfiles))
  nchem <- 1
  for (j in 1:length(allfiles)) {
    start <- regexpr("CAS_",allfiles[j])+4
    end   <- regexpr("_allstats.csv",allfiles[j])-1
    if(start>4 & end>start) {
      keep[j]<- TRUE
      allcasnums[nchem] <-substr(allfiles[j],start, end)
      nchem <- nchem+1
    }
  }
  keepfiles <- allfiles[keep]
  allfiledata <- list()
  for (j in 1:length(keepfiles)) {
    filedata<- read.csv(keepfiles[j])                                       #Load data from output file
    for (i in 1:length(metrics)) {
      datametric<-filedata[filedata$X==metrics[i],]                         # Identify metric of interest
      names(datametric)[2] <- "Cohort"                                      # Rename variables to add metric to name
      names(datametric)[3] <- paste("exp.dermal", metrics[i],sep="_")
      names(datametric)[4] <- paste("exp.ingest", metrics[i],sep="_")
      names(datametric)[5] <- paste("exp.inhal", metrics[i],sep="_")
      names(datametric)[6] <- paste("dose.inhal", metrics[i],sep="_")
      names(datametric)[7] <- paste("dose.intake", metrics[i],sep="_")
      names(datametric)[8] <- paste("abs.dermal.ug", metrics[i],sep="_")
      names(datametric)[9] <- paste("abs.ingest.ug", metrics[i],sep="_")
      names(datametric)[10] <- paste("abs.inhal.ug", metrics[i],sep="_")
      names(datametric)[11] <- paste("abs.tot.ug", metrics[i],sep="_")
      names(datametric)[12] <- paste("abs.tot.mgkg", metrics[i],sep="_")
      names(datametric)[13] <- paste("ddd.mass"  , metrics[i],sep="_")
      names(datametric)[14] <- paste("conc.max.prod.aer"  , metrics[i],sep="_")
      names(datametric)[15] <- paste("conc.max.prod.vap"  , metrics[i],sep="_")
      names(datametric)[16] <- paste("exp.inhal.indir"  , metrics[i],sep="_")
      if (i==1) {
        datametric <- datametric[,!names(datametric) %in% c("X")]
        CAS        <- allcasnums[j]
        alldata    <- cbind(CAS, datametric[ ,!names(datametric) %in% c("CAS")])   #Just reorder to put CAS in front
        setnames(alldata,names(alldata)[1],"CAS")
      } else {
        datametric <- datametric[,!names(datametric) %in% c("X", "Cohort")]        # Remove extraneous vars
        alldata    <- cbind(alldata,datametric)                                    # Add results for next metric
      }
      allfiledata[[j]]<-alldata
    }
    cat(paste("\n Processing chemical", j," of ",length(keepfiles)))
  }
  cat("\n Combining data... \n")
  finaldata <- (rbindlist(allfiledata))
  write.csv(finaldata,paste0("output/",run.name,"/",out.file), row.names=FALSE)
  return(finaldata)
}

#' filter_sources
#'
#' This function is used to select a subset of sources from another "source" file. 
#'
#' @param In.file a csv file that contains expsoure info default = "source_scen_prods.csv"
#' 
#' @param Out.file default = "source_scen_prods.csv"
#' 
#' @param IDS  default = "ALL"
#' 
#' @param Types default = "ALL"
#' 
#' @details This function is used to select a subset of sources from another "source" file.  This is achieved by specifying 
#' either a list of the desired source ids, or one or more sources types.  The allowed types are "A" for articles, "F" for
#' foods, or "P" for products.  Use the c() function to list more than one item (e.g. types=c("F","P") for foods and products).
#' Any of the three types of source files (that is, source_scen, source_chem, or source_vars) may be used.
#' Specifying specific sources is achieved using the "ids" argument.  This may require examining the in.file beforehand,
#' to obtain the correct source.id values for the desired sources.    
#' 
#' @return A .csv file of the same type as in.file, with the same variables and data, but fewer rows. The selected rows
#' have source.type matching one of the elements in the "types" argument, and also have source.id matching one of the
#' elements of the "ids" argument.  If either argument is missing, then all sources automatically match it.  Filter_sources
#' also returns an R object containing the same data as the output .csv file.
#' 
#' @export

filter_sources = function(in.file="source_scen_prods.csv",out.file="source_scen_subset.csv",ids=c("ALL"),types=c("ALL")){
  
  if (in.file==out.file)  {
    cat("\nInput and Output file have same name; must be different to avoid unintentionally overwrite.")
    stop()
  }
  sources <-read.csv(paste0("inputs/",in.file));
  setnames(sources,names(sources),tolower(names(sources)))
  if (types[1] !="ALL")  {
    types  <- trimws(tolower(types))
    stypes <- trimws(tolower(sources$source.type))
    sources<-sources[stypes %in% types,]
  }
  if (ids[1] !="ALL")  {
    ids    <- trimws(tolower(ids))
    sids   <- trimws(tolower(sources$source.id))
    sources<-sources[sids %in% ids,]
  }
  write.csv(sources, paste0("inputs/",out.file), row.names=FALSE)
  return(sources)
}




