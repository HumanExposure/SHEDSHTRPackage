#' A function used to export .rda files distributed with ShedsHT package
#'
#' This is a utility function used to export consumer product data sets into CSV format.
#'
#' @param data_set_name A string type parameter, represents name of an ".rda" data set. 
#'
#' @param output_pth A string type parameter, represents location to store exported CSV file. 
#'
#' @param output_fname A string type parameter, represents name of the exported CSV file.
#' 
#' @param quote_opt A string type parameter, represents a logical argument for if any character or factor 
#' columns should be surrounded by double quotes or not. 
#' 
#' @param na_opt A string type parameter, the string to use for missing values in the data.
#' 
#' @param row_names_opt A string type parameter, either a logical value indicating whether 
#' the row names of x are to be written along with x, or a character vector of row names to be written.
#' 
#' @return No variable will be returned. Instead, the function will save the data object into the file of choice. 
#' 
#' @import data.table
#' 
#' @import plyr
#' 
#' @import stringr
#' 
#' @import ggplot2
#' 
#' @export



ExportDataTables <- function(data_set_name, output_pth, output_fname, quote_opt=TRUE, na_opt = "", row_names_opt = FALSE){

	# first check if data set is available.
	if (!exists(data_set_name)){
		stop(paste0("Data set ", data_set_name, " is not available."))
	}

	output_pth_fname <- file.path(output_pth, output_fname)
	
	x <- eval(parse(text=data_set_name))
	if(tolower(substr(data_set_name,1,3))=="run") {
	  x <- as.data.table(x[,1])
	  setnames(x,names(x),"Variable,          Value")
	  quote_opt<-FALSE
	}  

	cc <- try(aa <- write.csv(x, output_pth_fname, quote=quote_opt, na=na_opt, row.names=row_names_opt), silent=T) 
	if(is(cc,"try-error")) {
		a <- "N/A"
		stop(paste0("Something is wrong! Please check output file path first."))
	} else{
	  	print(paste0(data_set_name, " is exported at ", output_pth_fname))
	}

}