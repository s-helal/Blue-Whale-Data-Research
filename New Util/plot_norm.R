
library(stringr)
library('zoo')
mval = 7

gold_standard_name <- 'C:/Users/samih/OneDrive - University of California, Davis/Percent Output Data/State/state_geo_3_tmp_1.csv'
goldstandard <- read.csv(gold_standard_name)
filepath <- 'C:/Users/samih/OneDrive - University of California, Davis/Percent Output Data/State'
setwd(filepath)

plot_norm <- function(xbad,xfixed,ybadfile,yfixedfile){
        ybad <- read.csv(ybadfile)[[1]]
        # yfixed <- read.csv(yfixedfile)[[1]]
        print(length(xbad))
        print(length(ybad))
        # make plot
        plot(xbad, rollmean(ybad,mval), ylim = c(0.2, 0.9),col="white", lwd=2, xlab = "", ylab = "", bty='l')
        lines(x4,rollmean(good,mval), col = "black", lwd=2) # 3km gold standard is always black
        lines(xbad,rollmean(ybad,mval), col = "red", lwd=2) # 8 km bad is red
        # lines(xfixed,rollmean(yfixed,mval), col = "darkcyan", lwd=2) # fixed is darkcyan
        
        # create legend
        ybadname <- str_remove(ybadfile, "_state")
        ybadname <- str_remove(ybadname, ".csv")
        ybadname <- str_replace(ybadname, '_', ' ')
        # yfixedname <- str_remove(yfixedfile, "_state")
        # yfixedname <- str_remove(yfixedname, ".csv")
        # yfixedname <- str_replace(yfixedname, '_', ' ')
        
        # legend("top", 
        #        legend = c("3 km 1 day", ybadname, yfixedname), 
        #        col = c("black", "red", "darkcyan"), 
        #        pch = c("-","-","-"), 
        #        bty = "n",
        #        pt.cex = 2,
        #        ncol = 3,
        #        cex =0.7, 
        #        text.width = 60,
        #        text.col = "black", 
        #        horiz = F , 
        #        xpd = TRUE,
        #        inset = c(-0.05, -0.3))
        
        # setup for savex
        dev.copy(png, str_replace(ybadfile,'_nofix.csv','.png'))
        dev.off()
}

x4 <- rollmean(seq(120,360, by=0.25)[1:959],mval)
x2 <- rollmean(seq(120,360, by=0.5)[1:479],mval)
x1 <- rollmean(seq(120,360, by=1)[1:239],mval)
good <- goldstandard[[1]]

# gold standard
bad <- 'state_geo_3_tmp_1.csv'
fixed <- 'state_geo_3_tmp_3_clouded_1.csv'
plot_norm(x4,x4,bad,fixed)

# 3 km 3 day
bad <- 'state_geo_3_tmp_3.csv'
fixed <- 'state_geo_3_tmp_3_clouded_1.csv'
plot_norm(x4,x4,bad,fixed)

# 3 km 8 day
bad <- 'state_geo_3_tmp_8.csv'
fixed <- 'state_geo_3_tmp_8_clouded_1.csv'
plot_norm(x4,x4,bad,fixed)

# 6 km 8 day
bad <- 'state_geo_6_tmp_8.csv'
fixed <- 'state_geo_6_tmp_8_clouded_1_rate_2.csv'
plot_norm(x4,x2,bad,fixed)

# 12 km 8 day (the worst case)
bad <- 'state_geo_12_tmp_8.csv'
fixed <- 'state_geo_12_tmp_8_clouded_1_rate_1.csv'
plot_norm(x4,x1,bad,fixed)
