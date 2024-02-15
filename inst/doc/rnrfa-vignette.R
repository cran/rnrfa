## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  # these should normally already be installed with rnrfa
#  install.packages(c("curl", "ggmap", "ggplot2", "httr", "jsonlite",
#                     "lubridate", "parallel", "sf", "tibble", "zoo"))

## ----eval=TRUE----------------------------------------------------------------
packs <- c("devtools", "DT", "leaflet", "dygraphs")
install.packages(packs, repos = "https://cloud.r-project.org")
lapply(packs, require, character.only = TRUE)

## ----eval=FALSE,message=FALSE-------------------------------------------------
#  install.packages("rnrfa")

## ----eval=FALSE---------------------------------------------------------------
#  devtools::install_github("ilapros/rnrfa")

## ----eval=FALSE,echo=TRUE-----------------------------------------------------
#  library(rnrfa)

## ----eval=TRUE,echo=FALSE-----------------------------------------------------
suppressPackageStartupMessages(library(rnrfa))

## -----------------------------------------------------------------------------
# Retrieve station identifiers:
allIDs <- station_ids()
head(allIDs)

## -----------------------------------------------------------------------------
# Retrieve information for all the stations in the catalogue:
allStations <- catalogue()
head(allStations)

## -----------------------------------------------------------------------------
# Define a bounding box:
bbox <- list(lon_min = -3.82, lon_max = -3.63, lat_min = 52.43, lat_max = 52.52)
# Filter stations based on bounding box
x <- catalogue(bbox)
dim(x); range(x$latitude); range(x$longitude)

# Filter based on minimum recording years
x <- catalogue(min_rec = 100)
dim(x); range(lubridate::year(x$`gdf-end-date`) - lubridate::year(x$`gdf-start-date`))

# Filter stations measuring a certain river
x <- catalogue(column_name="river", column_value="Wye")
dim(x); unique(x$river)

# Filter based on bounding box & metadata strings
x <- catalogue(bbox, column_name="river", column_value="Wye")
dim(x); unique(x$river)

# Filter stations based on threshold
x <- catalogue(bbox, column_name="catchment-area", column_value=">1")
dim(x); range(x$`catchment-area`)

# Filter based on minimum recording years
x <- catalogue(bbox, column_name = "catchment-area",
          column_value = ">1",
          min_rec = 30)
dim(x)

# Filter stations based on identification number
x <- catalogue(column_name="id", column_value="== c(3001,3002,3003)")
x$id

## -----------------------------------------------------------------------------
# Other combined filtering
someStations <- catalogue(bbox,
                          column_name = "id",
                          column_value = "==c(54022,54090,54091,54092,54097)",
                          min_rec = 35)

## -----------------------------------------------------------------------------
# Where is the first catchment located?
someStations$`grid-reference`$ngr[1]

# Convert OS Grid reference to BNG
osg_parse("SN853872")

## -----------------------------------------------------------------------------
# Convert BNG to WSGS84
osg_parse(grid_refs = "SN853872", coord_system = "WGS84")

## -----------------------------------------------------------------------------
osg_parse(grid_refs = someStations$`grid-reference`$ngr)

## ----fig.width=7,fig.asp=0.6--------------------------------------------------
# Fetch only time series data from the waterml2 service
info <- cmr(id = "3001")
plot(info)

## ----fig.width=7,fig.asp=0.6--------------------------------------------------
# Fetch time series data and metadata from the waterml2 service
info <- cmr(id = "3001", metadata = TRUE)
info$meta
plot(info$data, 
     main = paste("Monthly rainfall data for the", 
                info$meta$station.name,"catchment"), 
     xlab = "", ylab=info$meta$data.type.units)

## ----fig.width=7,fig.asp=0.6--------------------------------------------------
# Fetch only time series data
# info <- gdf(id = "3001")
# plot(info)

# Fetch time series data and metadata from the waterml2 service
info <- gdf(id = "3001", metadata = TRUE)
plot(info$data, 
     main = paste0("Daily flow data for the ", info$meta$station.name,
                   " catchment (",info$meta$data.type.units, ")"), 
     ylab = info$meta$data.type.name)

## ----fig.width=7--------------------------------------------------------------
# Search data/metadata
s <- cmr(c(3002,3003), metadata = TRUE)

# s is a list of 2 objects (one object for each site)
plot(s[[1]]$data, 
     main = paste(s[[1]]$meta$station.name, "and", s[[2]]$meta$station.name), 
     sub = "Catchment monthly rainfall", ylab = s[[1]]$meta$data.type.units)
lines(s[[2]]$data, col = "green")

s <- get_ts(c(3002, 3003), type = "gdf", metadata = TRUE)
plot(s[[1]]$data, 
     main = paste(s[[1]]$meta$station.name, "and", s[[2]]$meta$station.name), 
     sub = "Gauged daily flow", ylab=s[[1]]$meta$data.type.units)
lines(s[[2]]$data, col="pink2")

## -----------------------------------------------------------------------------
library(DT)
datatable(catalogue(column_name = "river", column_value = "Thames", all = FALSE))

## -----------------------------------------------------------------------------
library(leaflet)

leaflet(data = someStations) %>% addTiles() %>%
  addMarkers(~longitude, ~latitude, popup = ~as.character(paste(id,name)))

## -----------------------------------------------------------------------------
library(dygraphs)
dygraph(info$data) %>% dyRangeSelector()

## ----eval=FALSE---------------------------------------------------------------
#  library(parallel)
#  # Use detectCores() to find out many cores are available on your machine
#  cl <- makeCluster(getOption("cl.cores", detectCores()))
#  
#  # Filter all the stations within the above bounding box
#  someStations <- catalogue(bbox)
#  
#  # Get flow data with a sequential approach
#  system.time(s1 <- gdf(someStations$id, cl = NULL))
#  
#  # Get flow data with a concurrent approach (using `parLapply()`)
#  system.time(s2 <- gdf(id = someStations$id, cl = cl))
#  
#  stopCluster(cl)

## -----------------------------------------------------------------------------
# Linear model
library(ggplot2)
ggplot(allStations[!is.na(allStations$qmed),], 
       aes(x = as.numeric(`catchment-area`), y = qmed)) +
  geom_point() +
  stat_smooth(formula = y ~ x, method = "lm", col = "red") +
  xlab(expression(paste("Catchment area [Km^2]", sep=""))) +
  ylab(expression(paste("Mean flow [m^3/s]", sep="")))

