library(data.table)

custom_read = function(filename, header = T, colClasses = NULL, fill = F){
  if (grepl("gz$", filename)){
    return(fread(paste('zcat ', filename, sep=""), header = header, sep="\t", stringsAsFactors=F, showProgress=F, colClasses = colClasses, fill = fill))
  } else {
    return(fread(filename, header = header, sep="\t", stringsAsFactors=F, showProgress=F, colClasses = colClasses, fill = fill))
  }
}