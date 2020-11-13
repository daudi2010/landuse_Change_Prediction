# This script Outputs a csv of extracted values from the rasters
# Uses Brt to run  gbm functions to extract ROC and Dev  parameters using 25k points
# Computes relative influence of each brt iteration epoch
# plots pdp plots
# Revised May 2018
# Created by David Kanyari -daudi2010@gmail.com | Nancy Waithera
# All rights reserved

# set working #directory

setwd("E:\\NancyRStudyArea\\Final_datasets")
# Libraries
require(gbm) # load required libraries
require(dismo)
library(splancs) # For generating random points
library(foreign)
library(maptools)
library(spatstat)
library(dplyr)
require(sp)
library(rgdal)
library(raster)
library(csv)
library(rgeos)
#Additional libraries
source("brt.functions.R") # Brt functions script
source('utils.R')# file utilities..Important!
#

#Data folders
dependent <- c("dependent")# folder with dependent rasters
constants <-
  c("constants")# static data folder - Elevation, slope, roads, towns,rivers
#folders <- c("pov", "pop", "precip", "fires")# variable rasters  #original
folders <- c("pov2", "pop2", "precip", "fires")# variable rasters
#counties in the study area
counties<-c("study_area2.shp")
#--------------------------------------------#
ct<-file_path_sans_ext(counties[1] )# hardcode Study area
# county
c<-counties[1]
# create a #directory  for county
#dir.create(file.path("data", ct)) #  looping will be done later

print(counties[1])

# create a #directory for the log per county
#dir.create(file.path(paste("data", ct,sep="/"),"log"))

# craete data folder inside county file
#dir.create(file.path(paste("data", ct,sep="/"),"data"))
# log file
logFile<-paste(paste("data", ct,sep="/"),paste("log","Brt_log_file.txt",sep="/"),sep="/")
cat("BRT log file ... ", file=logFile, append=FALSE, sep = "\n")

countyf<-file.path("data", ct)
# select county from study area
county <- readOGR(dsn = paste("study_area",c,sep = "/"), layer = ct)

years <-
  
 c(
   #1995,
   #2000
  #,
   2005
   , 
   2010
   #2014
   ) # 1990 data not a complete set

cat ("Extracting values and creating csv for BRT")
for (i in 1:length(years)) {
  print(years[i])
  
  # This will be written after extraction of values
  extractedcsv <-
    paste(paste(countyf,"data",sep="/"),paste("/","Extractedvalues_"), paste(years[i], ".csv", sep = ""), sep = "")
  # check file exists
  
  if (!file.exists(extractedcsv)) {
    # Extract points
    cat("Reading existing points shapefile..\n")
   
    
    
    #Import SpatialPointsDataFrame
    pointspdf <-
      readOGR(dsn = "points/samplePoints.shp", layer = "samplePoints")
    
    #Intersect with study area poly# selected county  
    cat("intersecting points and polygons! study area and sample points")
    
    #pointspdf<-gIntersection(pointspdf1, county,
    #                         byid=FALSE, 
    #                         id=NULL,
    #                         drop_lower_td=FALSE, 
    #                         unaryUnion_if_byid_false=TRUE, 
    #                         checkValidity=FALSE)
    
    
    rastersC <- c() # a vector to hold  all rasters
    
    # add dependent variables per year
    print("Adding dependent variable")
    dependents <- list.files(path = dependent[1], pattern = "\\.tif$")
    dependents_img <-
      list.files(path = dependent[1], pattern = "\\.img$")
    # Loop and filter per year using grep..tiffs
    for (filed in grep(years[i], dependents, value = TRUE)) {
      #add to raster list build
      rastersC <-
        append(rastersC, paste(dependent[1], filed, sep = "/"))
      
    }
    #img
    for (filed in grep(years[i], dependents_img, value = TRUE)) {
      #add to raster list build
      rastersC <-
        append(rastersC, paste(dependent[1], filed, sep = "/"))
      
    }
    
    print("Adding constants independent variables")
    # add constants to rasters list build
    rasters <- list.files(path = constants[1], pattern = "\\.tif$")
    for (ras in rasters) {
      rastersC <- append(rastersC, paste(constants[1], ras, sep = "/"))
    }
    print("Adding other variable independnt variables")
    for (folder in folders) {
      # Extract only year needed
      for (fileg in list.files(path = folder,
                               pattern = paste(years[i], "\\.tif$", sep = ""))) {
        rastersC <- append(rastersC, paste(folder, fileg, sep = "/"))
        
      }
    }
    print (rastersC)
    # Extract points here!!
    # loop and extract
    cat('Extracting values to points..\n')
    data <- data.frame(coordinates(pointspdf))
    datanames <- c("x", "y")
    for (i in 1:length(rastersC)) {
      # Extract values to points....
      cat (paste("Extracting from", rastersC[i], sep = "  "))
      
      rastvals <- extract(
        raster(rastersC[i]),
        # pass a raster object!!!
        pointspdf,
        method = 'bilinear',
        small = FALSE,
        fun = mean,
        na.rm = TRUE
      )
      #column name.Create meaningful column names.Replace / with _
      col <- gsub("/", "_", rastersC[i])
      #col<- gsub("dependent", "", col)
      col <- gsub(".tif", "", col)
      col <- gsub(".img", "", col)
      col <- gsub("constants", "", col)
      print (col)
      # col <- file_path_sans_ext(file.names[i])
      
      #  output as dataframme
      cat('\nAdding  column to to Dataframe\n')
      #Add column  data to data frame
      data[[col]] <- rastvals
      # Add column  name to data frame
      datanames[2 + i] <- col
      cat("Done.. \n")
    }
    names(data) <- datanames
    
    
    cat('Sanitizing data..\n')
    # force {0.1} for dependent variables for consistency.A candidate bug in BRT
    
    data[, grepl("dependent", colnames(data))] <-
      round(data[, grepl("dependent", colnames(data))], 0)
    # remove outliers
    data[data == -9999] <- NA
    
    #preview data
    str(data)
    head(data)
    
    cat('Creating csv file\n')
    #output csv
    write.table(
      data,
      file = extractedcsv,
      sep = ",",
      col.names = NA,
      append = FALSE,
      qmethod = "double"
      
    )
    cat('Done..\n')
    
  } else{
    # end file exists
    
    #read existing file
    data <- read.csv(file = extractedcsv, as.is = T)
  }#end file exists block
  
  #define colums with dependent variables
  colums_ <- c(4, 5, 6)#df-mf,f-cl,f-gl  columns
 #-----------------------------------##
  # Run model on  25 iterations
  
  iteration_no=1
  
   # consider f-cl for trial
  #colums_ <- c(5)#df-mf,f-cl,f-gl  columns
  totalcolums <- ncol(data)# dependents start from colum 7--->
  
  for (col in 1:length(colums_)) {
    #df-mf,f-cl,f-gl  columns
    # select only those with value 1
    cat(paste(" \nExecuting ", colnames(data)[col]))
    
    
    cat("#######################################################################", file=logFile, append=TRUE, sep = "\n")
    
    
    
    colg <- colums_[col]
    # select only those with value 1
    gooddata <- data[data[[colg]] == 1, ]
    # get the no of rows
    rows_good <- nrow(gooddata)
    # remove  the  good rows from main data
    nulldata <- data[-rows_good, ]
    #  random sample of zeros  60:40 ratio
   sample <-
      sample.int(nrow(nulldata), floor(5.6 * nrow(gooddata) ), replace =
                  FALSE)
    # training data
    model.data <- nulldata[sample, ]
    
    # combine data frames
    
     new_mydata <- rbind(gooddata, model.data)
    # Randomize  the new dataframe rows
    model.data2 <- new_mydata[sample(nrow(new_mydata)),]
    #model.data2<-data  # portion of 100000 (noneed to split)
     # column  names for independent values only
    fd<-strsplit(names(data)[7:totalcolums]," ")
    #fd<-strsplit(names(data)[8:totalcolums]," ")
    #Vectors to hold results for  contributions
    
    # use all the 100000 data points!!!
    model.data2=data
    
    for ( f in 1:length(fd)){
      #global variables
      print(fd[[f]])
      
    }
    
    cat("done")
  
    
    
    filename<-paste(names(model.data2)[colg],".txt",sep="")
    
    
    cat(names(model.data2)[colg], file=logFile, append=TRUE, sep = "\n")
  
    mat <- matrix(, nrow = iteration_no+1, ncol = length(fd))
    for(i in 1:iteration_no){
      cat (paste(i," :iteration",sep = ""))
  
    response.tc5.lr005 <- gbm.step(
      data = model.data2,
      # data
      gbm.x = 8:totalcolums - 1, # fire 1995 has missing values
      #gbm.x = 8:totalcolums,
      #eliminate fires..too many missing values
      # colum 5 to 13
      gbm.y = colums_[col],
      # column 1
      family = "bernoulli",
      #define parameter for gbm.step function
      tree.complexity = 6,
      learning.rate = 0.005,
      bag.fraction = 0.75,
      step.size = 100  #  default  50
      #n.trees = 15000 # default 10000
    )
    # do relative influence
    
   
    
  for ( f in 1:length(names(data)[8:totalcolums-1])){
      #  get contribtions  and append values
      
       g<-6+f
      
       mat[i,f]<-response.tc5.lr005$contributions[names(data)[g],"rel.inf"]
       
       names(data)[g]
       response.tc5.lr005$contributions[names(data)[g],"rel.inf"]
       
      
    }
    }# end 25 loops
    
  
    # compute mean
    print("mean contribution")
    for ( f in 1:length(fd)){
      #global variables
      mat[iteration_no+1,f]<-mean(mat[,f],na.rm = TRUE)
      
      
    }
    print(colMeans(mat, na.rm = TRUE, dims = 1))
     # save csv file
    # row 26
    # create a ##directory for the plots per county
    ##dir.create(file.path(paste("data", ct,sep="/"),"contributions"))
    
    mat[iteration_no+1,]
    colnames(mat) <- fd
    mat
    write.table(mat,file=paste(paste("data", paste(ct,"contributions",sep="/"),sep="/"),filename,sep ="/"),append = FALSE) # keeps the rownames
    
    
    # create a #directory for the plots per county
    #dir.create(file.path(paste("data", ct,sep="/"),"plot"))
    # plot title
    title_<-paste(paste("data", ct,sep="/"),paste("plot","BRTplot_",sep="/"),sep="/")
    
    
     # save plots...Include column name in the plot
   # plotname <- paste(title_, colnames(data)[colg], sep = "_2_")
   # plotname <- paste(plotname, ".png", sep = "")
    
    
    # ROC
    roc2 <- mean(response.tc5.lr005$cv.roc.matrix)
    dev <-  response.tc5.lr005$cv.statistics$deviance.mean
    
    cat("mean deviance\n")
    print(dev)
    cat ("ROC\n")
    print(roc2)
    
    # output to log
    cat(paste("CV ROC",roc2,sep = "  : "), file=logFile, append=TRUE, sep = "\n")
    cat(paste("CV Deviance",dev,sep = "  : "), file=logFile, append=TRUE, sep = "\n")
    cat(paste("CV Deviance Se",response.tc5.lr005$cv.statistics$deviance.se,sep = "  : "), file=logFile, append=TRUE, sep = "\n")
    cat(paste("CV correlation",response.tc5.lr005$cv.statistics$correlation.mean,sep = "  : "), file=logFile, append=TRUE, sep = "\n")
    cat(paste("CV correlation Se",response.tc5.lr005$cv.statistics$correlation.se,sep = "  : "), file=logFile, append=TRUE, sep = "\n")
    #cat(paste("RelativeInf",summary(response.tc5.lr005),sep = "  : "), file=logFile, append=TRUE, sep = "\n")
    
    
    
    #png(filename=plotname)
    
    # plot pdp plots
   gbm.plot(response.tc5.lr005,
           n.plots = length(names(model.data2)[7:totalcolums-1]),
             # No ofplots vs independent variables
             write.title = FALSE)
    
    #dev.off()
    #windows()
    #savePlot(filename=plotname,device=dev.cur())
    # fitted values
    cat("Plot fitted values\n")
    
    #gbm.plot.fits(response.tc5.lr005)
    
    
    
    #cat("Interactions \n")
   
    # Relative influence
    #cat("Relative influence..\n")
    #cat(paste("Relative influence",summary(response.tc5.lr005),sep = "  : "), file=logFile, append=TRUE, sep = "\n")
    
    print(summary(response.tc5.lr005))
    
    
    
  }# end for loop - colums
  
  
}#  end epoch year loop!
