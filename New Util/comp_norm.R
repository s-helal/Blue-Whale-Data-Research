gold_standard_name <- 'C:/Users/samih/OneDrive - University of California, Davis/Percent Output Data/State/state_geo_3_tmp_1.csv'
goldstandard <- read.csv(gold_standard_name)

filepath <- 'C:/Users/samih/OneDrive - University of California, Davis/Percent Output Data/State'

library(stringr)
# set working directory
setwd(filepath)

patt = "*.csv"

file_names <- list.files(path=filepath, pattern=patt, full.names=FALSE, recursive=FALSE)
file_labels <- str_remove(file_names, "_yr_2008")
file_labels <- str_remove(file_labels, "state_output_data_")
file_labels <- str_remove(file_names, ".csv")

comp <- function(ser1, ser2){l2 <- sqrt(sum((ser1 - ser2)^2)); return(l2)}

indices <- 1:length(file_names)

for(i in indices){
        file_label <- file_labels[i]
        file <- read.csv(file_names[i])
        
        if(grepl('rate_1', file_label)){
                skip <- 4
                ind <- seq(1,dim(goldstandard)[1],by=skip)
                ind <- ind[0:-1]
        } else if(grepl('rate_2',file_label)){
                skip <- 2
                ind <- seq(1,dim(goldstandard)[1],by=skip)
                ind <- ind[0:-1]
        } else{
                skip <- 1
                ind <- seq(1,dim(goldstandard)[1],by=skip)
        }
        
        gold_standard <- goldstandard[ind,]
     
        file[is.na(file)] <- 0
        l2 <- comp(file, gold_standard)
        print(file_labels[i])
        # print(paste("rate", rate))
        print(l2)
        # print(dim(file)[1])
}

