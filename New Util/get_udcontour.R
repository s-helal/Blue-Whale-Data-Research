## load packages
require(adehabitatHR)
require(sp)
require(ggplot2)
require(plot3D)

## start timer
ptm = proc.time()
library(stringr)
# set working directory
filepath <-"C:/Users/samih/OneDrive - University of California, Davis/Percent Output Data/Sep"
setwd(filepath)
patt = "*.csv"

file_names <- list.files(path=filepath, pattern=patt, full.names=FALSE, recursive=FALSE)
file_labels <- str_remove(file_names, ".csv")
file_labels <- str_remove(file_labels, "sep_")

## function to load file into df with proper labels
prep_file <- function(file_name){
  # load file
  df_temp <- read.csv(file_name, header=FALSE)
  
  # add id column
  cols <- c("x","y","id")
  
  # format ids and colnames
  mylen <- dim(df_temp)[1]
  label <- rep(c(file_labels[i]), mylen)
  df_temp$id <- label
  colnames(df_temp) <- cols
  return(df_temp)
}

## load map
map_df <- read.csv("C:/Users/samih/OneDrive/Documents/Blue Whale Data Resolution Project/krill.csv", 
                   header = TRUE)
map_matrix <- t(as.matrix(map_df))

## set map to binary land/ocean
map_mask <- map_matrix
for(i in 1:length(map_mask)){
  map_mask[i] <- ifelse(is.na(map_mask[i]), 0, 1)
}

## prepare for loop
indices <- 1:length(file_names)
xvals <- seq(0, 359*3000, by = 3000)
yvals <- seq(0, 360*3000, by = 3000)
df_standard <- prep_file("sep_geo_3_tmp_1.csv")

## make sure the gold standard is first alphabetically
df_standard$id <- paste('0',df_standard$id[1], sep='')

## loop through files
for (i in indices){
  df_compare <- prep_file(file_names[i])
  
  # stack dataframes and rearrange columns
  df <- rbind(df_standard, df_compare)
  df <- df[,c("id","x","y")]
  
  # set up SpatialPointsDataFrame for the dataframe
  coordinates(df) <- c("x", "y")
  proj4string(df) <- CRS(as.character(NA))
  
  # compute utilization distributions and generate heatmap
  ud <- kernelUD(df, h = "href")
  
  # get ud contour
  ver_test <- getverticeshr(ud, percent = 60)
  
  # get x and y limits
  xmin <- min(df_compare$x, df_standard$x)
  xmax <- max(df_compare$x, df_standard$x)
  ymin <- min(df_compare$y, df_standard$y)
  ymax <- max(df_compare$y, df_standard$y)

  # set colors
  # red for gold standard, black for regular, mediumorchid for fixed, darkcyan for whales
  land_color = gray.colors(2, end = 1, alpha = .4)
  whale_color <- alpha('darkcyan', 0.3)
  # if(grepl('rate_4',df_compare$id[1])){
  #   ud_colors <- c('black', 'red2')}else{
  #   ud_colors <- c('black', 'mediumorchid')
  #   }
  ud_colors <- c('black', 'red2')
  
  # make plot
  png(paste(file_labels[i], "_ud.png", sep = ""))
  par(mar = c(2, 2, 2, 2))
  image2D(map_mask, yvals, xvals, col = land_color,
          ylab = NULL, main = NULL,
          axes = FALSE, frame.plot = FALSE, colkey = FALSE,
          title = file_labels[i])
  # points(df_compare$x, df_compare$y, col = whale_color, pch = 19)
  # plot(ver_test, lty = 1, lwd = c(3,3), ylim = c(ymin, ymax), xlim = c(xmin, xmax),
  #      border = ud_colors, add = TRUE)
  
  # allow images to show
  dev.off()
}

## end timer
proc.time() - ptm