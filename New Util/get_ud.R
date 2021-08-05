# start the clock
ptm <- proc.time()
library(stringr)
# set working directory
filepath <-"C:/Users/samih/OneDrive - University of California, Davis/Percent Output Data/Sep"
setwd(filepath)

# ## generate filenames for comparison
# geo_vals <- rep(paste(c(3,6,9,12)),3)
# tmp_vals <- rep(paste(c(1,3,8)),4)
# rate <- "4" # set to 4 if no rate
# 
# rate_string <- ifelse(rate == 4, "", paste("_rate_", rate, sep = ""))
# file_names <- paste("output_data_geo_", geo_vals, "_tmp_", tmp_vals, rate_string, "_sep", ".csv", sep = "")
# file_labels <- paste("geo_", geo_vals, "_tmp_", tmp_vals, rate_string, sep = "")

# patt = "^output_data_geo_\\d*_tmp_\\d_clouded_\\d_rate_\\d_sep.csv$"
patt = "*.csv"

file_names <- list.files(path=filepath, pattern=patt, full.names=FALSE, recursive=FALSE)
file_labels <- str_remove(file_names, ".csv")
file_labels <- str_remove(file_labels, "sep_")

# add gold standard
file_names <- c(file_names, "sep_geo_3_tmp_1.csv")
file_labels <- c(file_labels, "0gold_standard")

prep_file <- function(i){
  # load file
  df_temp <- read.csv(file_names[i], header=FALSE)
  
  # add id column
  cols <- c("x","y","id")
  
  # format
  mylen <- dim(df_temp)[1]
  label <- rep(c(file_labels[i]), mylen)
  df_temp$id <-  label
  colnames(df_temp) <- cols
  return(df_temp)
}

df <- prep_file(1)
## prepare for loop
indices <- 2:length(file_names)

## loop through files
for (i in indices){
  df_next <- prep_file(i)
  df <- rbind(df, df_next)
}


df <- df[,c("id","x","y")]

# set up SpatialPointsDataFrame for the dataframe
library(sp)
coordinates(df) <- c("x", "y")
proj4string(df) <- CRS(as.character(NA))

# compute utilization distributions and generate heatmap
library(adehabitatHR)
ud <- kernelUD(df, h = "href")

# image_2 <- image(ud)
udoverlap_2 <- kerneloverlap(df, method = "VI", percent = 95, conditional = FALSE)
# title(main = paste("steps/day =", rate, sep = " "), adj = 1, line = -15)

# for comparing to og, udoverlap_nov[,7]
main_uds <- udoverlap_2[,1]

for(item in main_uds){
  print(item)
}

# stop the clock
proc.time() - ptm