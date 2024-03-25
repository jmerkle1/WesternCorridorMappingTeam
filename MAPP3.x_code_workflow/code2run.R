
#packages
library(sf)
library(parallel)
library(raster)
library(terra)
library(devtools)

# load up functions from github
url <- "https://raw.githubusercontent.com/jmerkle1/WesternCorridorMappingTeam/main/MAPP3.x_code_workflow/functions"
funs <- c("CalcBBMM","CalcCTMM","CalcDBBMM","CalcKernel","CalcLineBuff",
          "CalcPopFootprint","CalcPopGrid","CalcPopUse","CalcSeqDistances")
for(i in funs){
  source_url(paste0(url,"/",i,".R"))
}
rm(url, funs)

# set WD to export folder from Migration Mapper 3.0
setwd("C:/Users/jmerkle_local/Desktop/Mapp_Practice/Moose_Jackson_App3_CMT/Exports")

# which seasons do you want to do?
dir("./sequences")  # there will be some values here. Put the ones you want in the next line
seasons <- c("Fall","Spring")



# ------------------#
# code for app 4 ####
# ------------------#

# ---------------------------- #
# 1. create overall raster ####
# ---------------------------- #
# Must do this based on ALL the data for the project, not just for one season !!!!

## First make d using all the formatted sequqnces for all the subherds
d <- readRDS("./workingFile.rds")
meta <- readRDS("./configOptions.rds")
d <- d[[2]]
meta$masterCrs
d <- st_as_sf(d, coords = c("x", "y"), 
              crs=meta$masterCrs)
head(d)

if(dir.exists("./PopulationGrid")==FALSE){
  dir.create("./PopulationGrid")
}

CalcPopGrid(datasf=d,
            out.fldr= "./PopulationGrid",
            mult4buff=0.3,
            cell.size=500)
rm(d)

# ------------------------- #
# 2. Calculate distances ####
# ------------------------- #

for(e in 1:length(seasons)){
  # read in the file
  d <- readRDS(paste0("./sequences/", seasons[e], "/", seasons[e], ".rds"))
  d <- st_transform(d, crs=meta$masterCrs)
          
  dists <- CalcSeqDistances(datasf=d, id.name="mig")
  head(dists)
    
  if(dir.exists("./Metadata")==FALSE){
    dir.create("./Metadata")
  }
          
  # write out:
  write.csv(dists, 
            file = paste0("./Metadata/", seasons[e], "_migration_distance_info.csv"),
            row.names=FALSE)
}


# ---------------------------------- #
# 3. Calculate UDs and Footprints ####
# ---------------------------------- #

Footprint.fldr <- "./Footprints"
UD.fldr <- "./UDs"

if(dir.exists(Footprint.fldr)==FALSE){
  dir.create(Footprint.fldr)
}
if(length(dir(Footprint.fldr))> 0)
  stop("Your Footprint.fldr Has something in it. It should be empty!")

if(dir.exists(UD.fldr)==FALSE){
  dir.create(UD.fldr)
}
if(length(dir(UD.fldr))> 0)
  stop("Your UD.fldr Has something in it. It should be empty!")

# choose 1 method for calculating UDs and/or footprints
## This needs to be a single method for the if() statements 
#opts <- c("LineBuff","BBMM", "DBBMM","Kernel","CTMM")
opts <- "BBMM"

# must first loop over seasons
for(e in 1:length(seasons)){
  
  # read in the file
  d <- readRDS(paste0("./sequences/", seasons[e], "/", seasons[e], ".rds"))
  d <- st_transform(d, crs=meta$masterCrs)
          
  # create the folders
  if(dir.exists(paste0(Footprint.fldr,"/",seasons[e]))==FALSE){
    dir.create(paste0(Footprint.fldr,"/",seasons[e]))
  }
  if(length(dir(paste0(Footprint.fldr,"/",seasons[e])))> 0)
    stop("Your paste0(Footprint.fldr,"/",seasons[e]) Has something in it. It should be empty!")
  
  if(dir.exists(paste0(UD.fldr,"/",seasons[e]))==FALSE){
    dir.create(paste0(UD.fldr,"/",seasons[e]))
  }
  if(length(dir(paste0(UD.fldr,"/",seasons[e])))> 0)
    stop("Your paste0(UD.fldr,"/",seasons[e]) Has something in it. It should be empty!")
  
  loopit <- unique(d$mig)  # you will loop over each id_yr_seas
  no_cores <- detectCores() - 5 # this should be the default, but the user could choose too
  
  # Setup cluster
  clust <- makeCluster(no_cores) 
  # export the objects you need for your calculations from your environment to each node's environment
  clusterExport(clust, varlist=c("d","Footprint.fldr","UD.fldr", "opts","loopit",
                                 "CalcBBMM","CalcLineBuff", "CalcDBBMM","CalcKernel",
                                 "CalcCTMM","seasons","e"))
  # i = 1
  result.tbl <- do.call(rbind, clusterApplyLB(clust, 1:length(loopit), function(i){
    
    
    # need library() here for the packages your calculations require for your calculations
    library(sf)
    library(terra)
    
    # grab the sequence of interest
    tmp <- d[d$mig==loopit[i],]
    
    # prep some params that carry across the functions
    mult4buff <- 0.3
    contour = 99
    max.timeout = 3600*12
    date.name = "date"
    Pop.grd = "./PopulationGrid/PopGrid_empty.tif"
    
    
    # which function to use?
    if(opts == "LineBuff"){
      return(CalcLineBuff(
        seq.sf = tmp,
        seq.name = loopit[i],
        date.name = date.name,
        Footprint.fldr = paste0(Footprint.fldr,"/",seasons[e]),
        Pop.grd = Pop.grd,
        buff = 300
      ))
    }
    
    if(opts == "BBMM"){
      return(CalcBBMM(
        seq.sf = tmp,
        seq.name = loopit[i],
        date.name = date.name,
        UD.fldr = paste0(UD.fldr,"/",seasons[e]),
        Footprint.fldr = paste0(Footprint.fldr,"/",seasons[e]),
        Pop.grd = Pop.grd,
        BMVar = NULL,
        location.error = 20,
        max.lag = 8,
        contour = contour,
        time.step = 5,
        mult4buff = mult4buff,
        max.timeout = max.timeout
      ))
    }
    
    if(opts == "DBBMM"){
      return(CalcDBBMM(
        seq.sf=tmp,
        seq.name=loopit[i],
        date.name=date.name,
        UD.fldr=paste0(UD.fldr,"/",seasons[e]),
        Footprint.fldr=paste0(Footprint.fldr,"/",seasons[e]),
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
    
    if(opts == "Kernel"){
      return(CalcKernel(
        seq.sf=tmp,
        seq.name=loopit[i],
        date.name=date.name,
        UD.fldr=paste0(UD.fldr,"/",seasons[e]),
        Footprint.fldr=paste0(Footprint.fldr,"/",seasons[e]),
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
        UD.fldr=paste0(UD.fldr,"/",seasons[e]),
        Footprint.fldr=paste0(Footprint.fldr,"/",seasons[e]),
        Pop.grd=Pop.grd,
        Information.Criteria="AIC",
        contour=contour,
        mult4buff=mult4buff,
        max.timeout=max.timeout
      ))
    }
    
    
  }))
  stopCluster(clust)   # you must stop the parallelization framework
  
  head(result.tbl)
  # View(result.tbl)
  # write out result.tbl
  write.csv(result.tbl, 
            file = paste0("./Metadata/", seasons[e], "_metadata.csv"),
            row.names = FALSE)
  
  print(paste0("Season ", seasons[e], " is done!"))
  
}  # end of loop over seasons

# have a quick look at the results
# load up all the Footprints
foots <- rast(dir(paste0(Footprint.fldr,"/",seasons[1]), full.names = TRUE))
plot(foots[[1]]) # plot one of them
plot(sum(foots)) # plot the sum of all of them

# load up all the UDs
UDs <- rast(dir(paste0(UD.fldr,"/",seasons[1]), full.names = TRUE))
plot(UDs[[1]])  # plot one of them
plot(sum(UDs))  # plot the sum of all of them

# -------------------#
# Population Outputs #
# -------------------#

# ------------------#
# code for app 5 ####
# ------------------#

if(dir.exists("./finalOutputs")==FALSE){
  dir.create("./finalOutputs")
}

# folder output name
out.fldr <- "./finalOutputs/Migration"

# create the needed folders and check them
if(dir.exists(out.fldr)==FALSE){
  dir.create(out.fldr)
}
if(length(dir(out.fldr))> 0)
  stop("Your out.fldr has something in it. It should be empty!")

if(dir.exists(paste0(out.fldr,"/","popUseMerged"))==FALSE){
  dir.create(paste0(out.fldr,"/","popUseMerged"))
}
if(length(dir(paste0(out.fldr,"/","popUseMerged")))> 0)
  stop("Your paste0(out.fldr,"/", popUseMerged) has something in it. It should be empty!")

if(dir.exists(paste0(out.fldr,"/","footPrintsMerged"))==FALSE){
  dir.create(paste0(out.fldr,"/","footPrintsMerged"))
}
if(length(dir(paste0(out.fldr,"/","footPrintsMerged")))> 0)
  stop("Your paste0(out.fldr,"/", footPrintsMerged) has something in it. It should be empty!")

# calculate population use final products
CalcPopUse(
  UD.fldr = UD.fldr, 
  out.fldr = paste0(out.fldr, "/popUseMerged"), 
  seas2merge = seasons, 
  udFootprintsToDrop = NULL,
  merge.order = c("id","year"),  
  contour = 99,   
  contour.type = "Area",    
  contour.levels = c(5,10,15,20,30,40,50,60,70,80,90),  
  min_area_drop = 20000,  
  min_area_fill = 20000,  
  simplify = TRUE,     
  ksmooth_smoothness = 2,  
  out.proj = meta$masterCrs  
)

# calculate population footprint final products
CalcPopFootprint(
  Foot.fldr = Footprint.fldr, 
  out.fldr = paste0(out.fldr, "/footPrintsMerged"), 
  udFootprintsToDrop = NULL,
  seas2merge = seasons, 
  contour.levels = c(5,10,15,20,30),  
  min_area_drop = 20000,  
  min_area_fill = 20000,  
  simplify = TRUE,     
  ksmooth_smoothness = 2,  
  out.proj = meta$masterCrs  
)

# Check your folders! THey should match the folders that Migration Mapper 3.0 would do.

# have a look at what you did
Foot.conts <- st_read(paste0(out.fldr, "/footPrintsMerged"),
                      "Footprint_contours")
head(Foot.conts)
Foot.conts <- Foot.conts[order(Foot.conts$contour),]
table(Foot.conts$contour)
plot(Foot.conts)

Use.conts <- st_read(paste0(out.fldr, "/popUseMerged"),
                      "Pop_use_contours")
head(Use.conts)
table(Use.conts$contour)
plot(Use.conts)
