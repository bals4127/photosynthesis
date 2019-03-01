#' @title Reads and cleans data from PICARRO output file.
#'
#' @description Reads and cleans data from PICARRO output file.
#' @param file_path file path for picarro .dat file.
#' @return Reads and cleans data from PICARRO output file.
#' @importFrom lubridate with_tz
#' @importFrom lubridate parse_date_time
#' @importFrom utils read.delim
#' @export
read.picarro <- function (file_path) {
    picarro_in <- read.delim(file_path , sep = "",header = T)
    options(digits.secs = 3)
    picarro_in$DATE <- format(picarro_in$DATE, format = "%m/%d/%Y" , origin = "GMT")
    picarro_in$TIME <- format(picarro_in$TIME, format = "%H:%M:%S")
    original.timestamp=format(as.POSIXct(paste(picarro_in$DATE, picarro_in$TIME)), "%Y/%m/%d %H:%M:%OS")
    picarro.time = parse_date_time(original.timestamp, "%Y/%m/%d %H:%M:%OS" )
    picarro_to_pullman = with_tz(picarro.time,tzone = Sys.timezone())
    picarro_in$LocalTime = picarro_to_pullman
    return(picarro_in)
    }
