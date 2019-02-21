#' @title Reads and cleans data from TDL output file (.dat).
#'
#' @description Read, clean and compute C and O isotopic discrimination from TDL output file.
#' @param file_path file path for TDl .dat file.
#' @return Read, clean and compute C and O isotopic discrimination from TDL output file.
#' @export
#' @importFrom utils read.csv
read.tdl<- function(file_path){
  mydata <- read.csv(file_path, skip=4, header=F)
  h <- readLines(file_path)[1:4]
  tdl_colnames <- gsub("\"","",strsplit(paste(h, collapse=""), ",")[[1]])
  tsloc <- grep("timestamp",tolower(tdl_colnames))
  list_colnames <- tdl_colnames[tsloc:(tsloc+ncol(mydata)-1)]
  list_colnames<-gsub("[(,)]","",list_colnames)
  names(mydata) <- list_colnames
  names(mydata)[1] <- "DateTime"
  mydata$DateTime <- as.POSIXct(mydata$DateTime, tz = Sys.timezone()) #This will work if the TDL and your computer are in the same time zone.
  return(mydata)
}
