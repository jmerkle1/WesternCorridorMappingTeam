
library(sf)
library(parallel)
library(raster)

setwd("C:/Users/jmerkle_local/Desktop/new_mapp")

# create a practice dataset
# library(foreign)
# proj_of_dbfs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
# seqs_fldr <- "C:/Users/jmerkle_local/Desktop/Pronghorn_Sublette_phaseII/sequences_SE"
# fls <- list.files(seqs_fldr, ".dbf$")
# d <- do.call(rbind, lapply(1:length(fls), function(i){
#   db <- read.dbf(paste(seqs_fldr, fls[i],sep="/"), as.is=TRUE)
#   db$mig <- sub(".dbf","",fls[i])
#   return(db)
# }))
# d$date <- as.POSIXct(strptime(d$date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")  
# d <- st_as_sf(d, coords = c("x","y"), dim="XY", crs=proj_of_dbfs)
# saveRDS(d, paste0(getwd(), "/data4practice.rds"))

# load up a practice dataset
d <- readRDS("data4practice.rds")

# ------------------#
# code for app 4 ####
# ------------------#


# ---------------------------- #
# 1. create overall raster ####
# ---------------------------- #
# Must do this based on ALL the data for the project, not just for one season !!!!
# Also must remove the mortality and bad points!

source("./functions/CalcPopGrid.R")
CalcPopGrid(datasf=d,
            out.fldr=getwd(),
            mult4buff=0.3,
            cell.size=500)
rm(CalcPopGrid)

# ------------------------- #
# 2. Calculate distances ####
# ------------------------- #
source("./functions/CalcSeqDistances.R")
dists <- CalcSeqDistances(datasf=d, id.name="mig")
head(dists)
# could write this out:
# write.csv(mig_dists, file=paste(metadata_fldr,"/migration_distance_info.csv",sep=""), row.names=FALSE)
rm(CalcSeqDistances)

# ---------------------------------- #
# 3. Calculate UDs and Footprints ####
# ---------------------------------- #
source("./functions/CalcBBMM.R")
source("./functions/CalcDBBMM.R")
source("./functions/CalcKernel.R")
source("./functions/CalcLineBuff.R")
source("./functions/CalcCTMM.R")

UD.fldr="./UDs"
Footprint.fldr="./Footprints"

#check the new directories
if(dir.exists(UD.fldr)==FALSE){
  dir.create(UD.fldr)
}
if(length(dir(UD.fldr))> 0)
  stop("Your UD.fldr Has something in it. It should be empty!")
if(dir.exists(Footprint.fldr)==FALSE){
  dir.create(Footprint.fldr)
}
if(length(dir(Footprint.fldr))> 0)
  stop("Your Footprint.fldr Has something in it. It should be empty!")

# method for calculating UDs and/or footprints
opts <- c("LineBuff","BBMM","dBBMM","kernel","CTMM")
opts <- "LineBuff"

loopit <- unique(d$mig)  # you will loop over each id_yr_seas
no_cores <- detectCores() - 1 # this should be the default, but the user could choose too

# Setup cluster
clust <- makeCluster(no_cores) 
# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("d","Footprint.fldr","opts","UD.fldr","loopit",
                               "CalcBBMM","CalcDBBMM","CalcKernel","CalcLineBuff","CalcCTMM"))

result.tbl <- do.call(rbind, clusterApplyLB(clust, 1:length(loopit), function(i){
  
  
  # need library() here for the packages your calculations require for your calculations
  library(sf)
  library(parallel)
  library(raster)
  
  # grab the sequence of interest
  tmp <- d[d$mig==loopit[i],]
  
  # prep some params that carry across the functions
  mult4buff <- 0.3
  Pop.grd="./PopGrid_empty.tif"
  contour=99
  max.timeout=3600*12
  date.name="date"
  UD.fldr=UD.fldr
  Footprint.fldr=Footprint.fldr
  
  # use new functions!!!!
  if(opts == "LineBuff"){
    return(CalcLineBuff(
      seq.sf=tmp,
      seq.name=loopit[i],
      date.name=date.name,
      Footprint.fldr=Footprint.fldr,
      Pop.grd=Pop.grd,
      buff=200
    ))
  }
  
  if(opts == "BBMM"){
    return(CalcBBMM(
      seq.sf=tmp,
      seq.name=loopit[i],
      date.name=date.name,
      UD.fldr=UD.fldr,
      Footprint.fldr=Footprint.fldr,
      Pop.grd=Pop.grd,
      BMVar=NULL,
      location.error=20,
      max.lag=8,
      contour=contour,
      time.step=5,
      mult4buff=mult4buff,
      max.timeout=max.timeout
    ))
  }
  
  if(opts == "dBBMM"){
   return(CalcDBBMM(
     seq.sf=tmp,
     seq.name=loopit[i],
     date.name=date.name,
     UD.fldr=UD.fldr,
     Footprint.fldr=Footprint.fldr,
     Pop.grd=Pop.grd,
     location.error=20,
     max.lag=8,
     contour=contour,
     dbbmm.margin=11,
     dbbmm.window=31,
     mult4buff=mult4buff,
     max.timeout=max.timeout
   )) 
  }

  if(opts == "kernel"){
    return(CalcKernel(
      seq.sf=tmp,
      seq.name=loopit[i],
      date.name=date.name,
      UD.fldr=UD.fldr,
      Footprint.fldr=Footprint.fldr,
      Pop.grd=Pop.grd,
      smooth.param=NULL,
      contour=contour,
      mult4buff=mult4buff,
      subsample=NULL,
      max.timeout=max.timeout
    ))
  }
  
  if(opts == "CTMM"){
    return(CalcCTMM(
      seq.sf=tmp,
      seq.name=loopit[i],
      date.name=date.name,
      UD.fldr=UD.fldr,
      Footprint.fldr=Footprint.fldr,
      Pop.grd=Pop.grd,
      Information.Criteria="AIC",
      contour=contour,
      mult4buff=mult4buff,
      max.timeout=max.timeout
    ))
  }

}))
stopCluster(clust)   # you must stop the parallelization framework


View(result.tbl)
# write out result.tbl
write.csv(result.tbl, file = "metadata.csv", row.names = FALSE)

# have a look at the results
foots <- stack(dir(Footprint.fldr, full.names = TRUE))
raster::plot(sum(foots))
# plot(d$geometry, add=T, pch=".")

UDs <- stack(dir(UD.fldr, full.names = TRUE))
raster::plot(mean(UDs))
# plot(d$geometry, add=T, pch=".")

# ------------------#
# code for app 5 ####
# ------------------#

source("./CalcPopUse.R")
source("./CalcPopFootprint.R")

CalcPopUse(
  UD.fldr = "C:/Users/jmerkle_local/Desktop/Corridor_test/UDs_Test", 
  out.fldr = "C:/Users/jmerkle_local/Desktop/Corridor_test/UD_out", 
  seas2merge = c("spring", "fall"), 
  udFootprintsToDrop = NULL,
  merge.order = c("year", "id"),  
  contour = 99,   
  contour.type = "Area",    
  contour.levels = c(5,10,15,20,30,40,50,60,70,80,90),  
  min_area_drop = 20000,  
  min_area_fill = 20000,  
  simplify = TRUE,     
  ksmooth_smoothness = 2,  
  out.proj = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"  
)

CalcPopFootprint(
  Foot.fldr = "C:/Users/jmerkle_local/Desktop/Corridor_test/Footprints_Test", 
  out.fldr = "C:/Users/jmerkle_local/Desktop/Corridor_test/Footprints_out", 
  udFootprintsToDrop = NULL,
  seas2merge = c("spring", "fall"), 
  contour.levels = c(5,10,15,20,30),  
  min_area_drop = 20000,  
  min_area_fill = 20000,  
  simplify = TRUE,     
  ksmooth_smoothness = 2,  
  out.proj = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"  
)


