##################################################################################
### Overview: This script takes all files in a folder and combines them into one
###            single file called combined data.csv with the filename as a row name
###
### author:   Alessandro Donada (alessandro.donada@gmail.com)
### date:     09/12/2022
###
##################################################################################


#set the working directory where all of your analysis files are
setwd("C:/Users/Desktop/Exported files")

#get the names of all the analysis files in the directory
filenames <- list.files(path = ".", pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)

#using the plyr library read in all of the individual files and combine them together using the rbind function which combines vectors by rows
library(plyr)
import.list <- llply(filenames, read.csv)
combined <- do.call("rbind", sapply(filenames, read.csv, simplify = FALSE))

#once the data have been merged into a single data frmae save it to a .csv file
write.csv(combined, file="combined data.csv")



  
