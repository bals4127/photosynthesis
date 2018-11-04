#' @title Read LI6400 file from csv format.
#'
#' @description Reads Li6400 file from csv format, removes unnessory coulmn, adds Li6400 name, measurment data and remarks in separate column.
#' @param file_path path of Li6400 csv file.
#' @return Read and clean csv output of Li6400 system. Adds Li6400 name, measurment date and Remarks (with respective observation log) as an extra column
#' @export
#' @importFrom lubridate mdy_hms
#' @importFrom utils read.csv
read.li6400<- function (file_path){
  ##Read remark function
  read_remark<- function(read2){
  rmk_df<- read2[read2$Obs == "Remark=", ]## remove remarks
  temp_obs<- as.numeric(as.character(droplevels
  (read2[as.numeric(rownames(rmk_df))+1, "Obs"])))
  remark<- cbind(temp_obs, rmk_df)[ , c("temp_obs", "HHMMSS")]
  names(remark)<- c("Obs","Remark" )
  return(remark)
  }
read1 <- read.csv(file_path)
where_do_i_start <- which(read1[,1]=="in") +1
# The number of lines before the data varies with each licor.
the_end <- nrow(read1) ## where I should end
open_ver<- colnames(read1)[1]
aa<-data.frame(droplevels(read1[[open_ver]][1]))
# Adding a timezone (computer not Li measurments locationfor clarity
datetime<- mdy_hms(aa[1,],tz=Sys.timezone())
licor.date<- as.Date(datetime)
start_time<- format(datetime, format= "%H:%M:%S")
#as.POSIXct(paste(licor.date, start_time), format="%Y-%m-%d %H:%M:%S")
##Licor name
licor_name<- levels(droplevels(read1[[2]][2]))
read2  <- read.csv(file_path, skip = (where_do_i_start - 2), header = T)
##remarks
remark<- read_remark(read2)
#read2<- merge(read2, remark, by="Obs", all.x=T)
read2  <- read2[!read2$Obs == "Remark=", ]## remove remarks
read2 = read2 [-1,]## bad row in data
for(i in c(5:ncol(read2))) {
read2[,i] <- as.numeric(as.character(read2[,i]))} ## remove the two-line header from Licor output, re-define columns as numeric
read2  <- droplevels(read2)
read2$datetime <- as.POSIXct(paste(licor.date,read2$HHMMSS, Sys.timezone()))
read2$licor<- licor_name
read2$Date<- licor.date
read2$Start_time<- start_time
## add remorks
test_r<- merge(read2, remark, by="Obs", all.x=T)
out<- test_r[order(test_r$datetime),]
return(out)
}
