#'allstats.variable.rank.plot
#'
#'A ggplot scatter plot with output.variable on the y axis and chemical rank on the x axis.
#'User can color points by cohort (see cohort.col argument) and assign percentiles to be points (see metrics argument)
#'
#'
#'@param run.name default = name of last run (in current R session)
#'
#'@param output.variable default = abs.tot.mgkg (mg/kg/day).
#'
#'@param metrics argument is a vector of character strings with desired percentiles or metrics.
#'
#'@param cohort.col argument vector of different cohorts. default = "Total". To add cohorts use a vector of concatenated character strings, such as cohort.col = c("Total", "Male", "Female")
#'
#'@return a plot of output.variable vs chemical rank
#'
#'@details A SHEDS run creates an output file for each chemical. This plot pulls from the allstats.csv for each chemical in the output file.
#'Chemicals are ranked in ascending order of output.variable on the x axis. The y axis is log transformed and a constant of 1e-10 added to y variable for
#'log transformation of 0's. The default color is black for the "Total" cohort, but user can include as many cohorts as they wish as long as they are calculated in the SHEDS run.
#'Percentiles must be calculated in SHEDS run.the selected cohorts.
#'
#'@export
#'

allstats.variable.rank.plot<- function(run.name=specs$run.name, output.variable = "abs.tot.mgkg" , metrics = c("50%", "95%", "mean"), cohort.col = "Total", label_chem = F) {
  allfiles   <-list.files(path=paste0("output/",run.name), all.files=FALSE, include.dirs=FALSE, full.names=TRUE)

   if(length(allfiles)<1){
    stop(paste0(run.name, " is not in the output folder."))
  }

  alldtxsid <- list()
  keep  <- rep(FALSE,length(allfiles))
  nchem <- 1
  for (j in 1:length(allfiles)) {
    start <- regexpr("CAS",allfiles[j])
    end   <- regexpr("_allstats.csv",allfiles[j])-1
    if(start>4 & end>start) {
      keep[j]<- TRUE
      alldtxsid[nchem] <-substr(allfiles[j],start, end)
      nchem <- nchem+1
    }
  }
  keepfiles <- allfiles[keep]
  allfiledata <- list()
  for (j in 1:length(keepfiles)) {
    filedata<- read.csv(keepfiles[j])
    pracfile2 <- list()
    for(i in 1:length(metrics)){
      datametric<-filedata[filedata$X == metrics[i],]
          if(nrow(datametric)<1){
            stop(paste0("The metric ", metrics[i]," was not calculated by the SHEDS-HT run."))
          }
      names(datametric)[1] <- "Statistic"
      names(datametric)[2] <- "Cohort"
      pracfile2 <- rbind(pracfile2,datametric)
      DTXSID        <- alldtxsid[j]
      alldata2    <- cbind(DTXSID, pracfile2)
      setnames(alldata2,names(alldata2)[1],"DTXSID")
    }
    cat(paste0("\n processing chemical ", j, " of ", length(keepfiles)))
    allfiledata[[j]]<-alldata2

  }
  cat("\n Combining data... \n")
  finaldata <- (rbindlist(allfiledata))
  finaldata <-as.data.frame(finaldata)

  finaldata_sub <- finaldata[,c("DTXSID", "Statistic", "Cohort", output.variable)]
  colnames(finaldata_sub) <- c("DTXSID", "Statistic", "Cohort", "output")

  ##Begin Visualization

  #What will the label on the y-axis be?
  if(output.variable == "exp.dermal")       {lab = "exp.dermal (ug/day)"}
  if(output.variable == "exp.ingest")       {lab = "exp.ingest (ug/day)"}
  if(output.variable == "exp.inhal")        {lab = "exp.inhal (ug/m3)"}
  if(output.variable == "dose.inhal")       {lab = "dose.inhal (ug/day)"}
  if(output.variable == "dose.intake")      {lab = "dose.intake (mg/kg/day)"}
  if(output.variable == "abs.dermal.ug")    {lab = "abs.dermal.ug (ug/day)"}
  if(output.variable == "abs.ingest.ug")    {lab = "abs.ingest.ug (ug/day)"}
  if(output.variable == "abs.inhal.ug")     {lab = "abs.inhal.ug (ug/day)"}
  if(output.variable == "abs.tot.ug")       {lab = "abs.tot.ug (ug/day)"}
  if(output.variable == "abs.tot.mgkg")     {lab = "abs.tot.mgkg (mg/kg/day)"}
  if(output.variable == "ddd.mass")         {lab = "ddd.mass (g.day)"}
  if(output.variable == "conc.max.prod.aer"){lab = "conc.max.prod.aer (ug/m3)"}
  if(output.variable == "conc.max.prod.vap"){lab = "conc.max.prod.vap (ug/m3)"}
  if(output.variable == "exp.inhal.indir")  {lab = "exp.inhal.indir (ug/m3)"}
  if(output.variable == "abs.hm.ug")        {lab = "abs.hm.ug (ug/day)"}



  #Create Graphs
  theme_set(theme_bw())

  ##Outdated code when there was a NULL option. Not Total is automatically the default
  if(is.null(cohort.col)){
    fc_null <-finaldata_sub[finaldata_sub$Cohort == "Total",]
    all_plot<-ggplot(fc_null, aes(x=as.numeric(reorder(DTXSID,as.numeric(output)))))+
      scale_y_log10()+
      labs(y = lab, x = "Chemical Rank")+
      geom_point(aes(y=as.numeric(output)+1e-10, shape = Statistic), col="black")+
      scale_shape_discrete(solid=F) + theme_bw()
    cat("\n Creating graph for all cohorts...\n")
    print(all_plot)
  }

  if(!is.null(cohort.col)){
    color_data <- data.frame()
    for(i in 1:length(cohort.col)){
    fc <- finaldata_sub[finaldata_sub$Cohort == cohort.col[i],]
    color_data<-rbind(color_data,fc)}
    
    color_data$x<-as.numeric(reorder(color_data$DTXSID,as.numeric(color_data$output)))
    color_data$y<-as.numeric(color_data$output)+1e-10
  
  cbpalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
  color_plot<-ggplot(color_data, aes(x=x, y=y,col = Cohort, shape = Statistic))+
    geom_point()+
    scale_y_log10()+
    labs(y = lab, x = "Chemical Rank")+
    scale_shape_discrete(solid=F)+
    scale_colour_manual(values = cbpalette)
  
  if (label_chem==T){
    
    test2<-color_data[which(color_data$Cohort=="Total" & color_data$Statistic=="95%"),]
    color_plot<-color_plot+geom_text(data=test2, aes(x=x,y=y,label=DTXSID), vjust=1, color="grey", angle=90,show.legend = FALSE)
  }   
    
}
  cat("\n Creating graph for specific cohorts...\n")
  print(color_plot)
}

#' puc.rank.plot is a ggplot scatter plot with desired srcMeans variable on the y axis and chemical rank on the x axis.
#' The dependent variable on the y axis can be changed using output.variable argument, but only one output variable can be plotted at a time.
#' User can either plot the individual product or the higher level puc using the combine argument. Data comes from all_srcMeans output file, so data points are the mean for that output variable.
#'
#'@param run.name default = name of last run (in current R session)
#'
#'@param output.variable argument is a single character string. default = "exp.dermal".
#'
#'@return a plot of output.variable vs chemical rank with a legend for color and shape
#'
#'@details A SHEDS run creates an output file for each chemical. This plot pulls from the srcMeans.csv for each chemical in the output file.
#'Chemicals are ranked in ascending order of output.variable on the x axis. The y axis is log transformed and a constant of 1e-10 added to y variable for
#'log transformation of 0's. The plot colors points by higher level PUC.
#'
#'@export

puc.rank.plot <- function(run.name, output.variable = "exp.dermal", label_chem = FALSE){
  allfiles   <-list.files(path=paste0("output/",run.name), all.files=FALSE, include.dirs=FALSE, full.names=TRUE)

  if(length(allfiles)<1){
    stop(paste0(run.name, " is not in the output folder."))
  }

  alldtxsid <- list()
  keep  <- rep(FALSE,length(allfiles))
  nchem <- 1
  for (j in 1:length(allfiles)) {
    start <- regexpr("CAS",allfiles[j])
    end   <- regexpr("_all_srcMeans.csv",allfiles[j])-1
    if(start>4 & end>start) {
      keep[j]<- TRUE
      alldtxsid[nchem] <-substr(allfiles[j],start, end)
      nchem <- nchem+1
    }
  }
  keepfiles <- allfiles[keep]
  allfiledata <- list()
  for (j in 1:length(keepfiles)) {
    filedata<- read.csv(keepfiles[j])
    filedata <- filedata[,-1] #remove first row -
    DTXSID        <- alldtxsid[j]
    filedata2    <- cbind(DTXSID, filedata)
    setnames(filedata2,names(filedata2)[1],"DTXSID")
    allfiledata[[j]]<-filedata2
    cat(paste0("\n processing chemical ", j, " of ", length(keepfiles)))
  }
  cat("\n Combining data... \n")
  finaldata <- (rbindlist(allfiledata))
  finaldata <- as.data.frame(finaldata)


  finaldata <- finaldata[,c("DTXSID", "src.names",  output.variable)]
  colnames(finaldata) <- c("DTXSID", "src.names", "output")

  #Establish Plot labels
  if(output.variable == "exp.dermal")       {lab = "exp.dermal (ug/day)"}
  if(output.variable == "exp.ingest")       {lab = "exp.ingest (ug/day)"}
  if(output.variable == "exp.inhal")        {lab = "exp.inhal (ug/m3)"}
  if(output.variable == "dose.inhal")       {lab = "dose.inhal (ug/day)"}
  if(output.variable == "f.dermal")         {lab = "f.dermal"}
  if(output.variable == "f.ingest")         {lab = "f.ingest"}
  if(output.variable == "f.inhal")          {lab = "f.inhal"}
  if(output.variable == "mean.mass")        {lab = "mean.mass (ug)"}

  options(dplyr.summarise.inform = FALSE)

    plotdata <- finaldata %>%
      mutate(Category = ifelse(src.names == "Total", "Total" , src.names))

  theme_set(theme_bw())
  cbpalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999", "#E411E0")


    #Create graph of output.variable level by PUC
  cat("\n Creating puc plot...\n")
  
  plotdata$x<-as.numeric(reorder(plotdata$DTXSID,as.numeric(plotdata$output)))
  plotdata$y<-as.numeric(plotdata$output)+1e-10
  
  puc.path.plot <- ggplot(plotdata, aes(x=x,y=y,color = Category))+
      geom_jitter(width=0.2)+
      scale_y_log10()+
      #scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
      labs(y = paste0("Mean ", lab), x = "Chemical Rank", color = "PUC")+
      scale_colour_manual(values = cbpalette)
  
      if (label_chem==T){
        
        test2<-plotdata[which(plotdata$Category=="Total"),]
        puc.path.plot<-puc.path.plot+geom_text(data=test2, aes(x=x,y=y,label=DTXSID), vjust=1, color="grey", angle=90,show.legend = FALSE)
      }   
  
  
  print(puc.path.plot)

}


