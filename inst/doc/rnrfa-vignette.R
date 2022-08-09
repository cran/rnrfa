## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)

## -----------------------------------------------------------------------------
#  install.packages(c("cowplot", "httr", "xts", "ggmap", "ggplot2", "sp", "rgdal", "parallel", "tibble"))

## -----------------------------------------------------------------------------
#  packs <- c("devtools", "DT", "leaflet")
#  install.packages(packs)
#  lapply(packs, require, character.only = TRUE)

## -----------------------------------------------------------------------------
#  install.packages("rnrfa")

## -----------------------------------------------------------------------------
#  devtools::install_github("ilapros/rnrfa")

## -----------------------------------------------------------------------------
#  library(rnrfa)

## -----------------------------------------------------------------------------
#  # Retrieve station identifiers:
#  allIDs <- station_ids()
#  head(allIDs)

## -----------------------------------------------------------------------------
#  # Retrieve information for all the stations in the catalogue:
#  allStations <- catalogue()
#  head(allStations)

## -----------------------------------------------------------------------------
#  # Define a bounding box:
#  bbox <- list(lon_min = -3.82, lon_max = -3.63, lat_min = 52.43, lat_max = 52.52)
#  # Filter stations based on bounding box
#  catalogue(bbox)

## -----------------------------------------------------------------------------
#  # Filter based on minimum recording years
#  catalogue(min_rec = 100)
#  
#  # Filter stations belonging to a certain hydrometric area
#  catalogue(column_name="river", column_value="Wye")
#  
#  # Filter based on bounding box & metadata strings
#  catalogue(bbox, column_name="river", column_value="Wye")
#  
#  # Filter stations based on threshold
#  catalogue(bbox, column_name="catchment-area", column_value=">1")
#  
#  # Filter based on minimum recording years
#  catalogue(bbox, column_name = "catchment-area",
#            column_value = ">1",
#            min_rec = 30)
#  
#  # Filter stations based on identification number
#  catalogue(column_name="id", column_value=c(3001,3002,3003))

## -----------------------------------------------------------------------------
#  # Other combined filtering
#  someStations <- catalogue(bbox,
#                            column_name = "id",
#                            column_value = c(54022,54090,54091,54092,54097),
#                            min_rec = 35)

## -----------------------------------------------------------------------------
#  # Where is the first catchment located?
#  someStations$`grid-reference`$ngr[1]
#  
#  # Convert OS Grid reference to BNG
#  osg_parse("SN853872")

## -----------------------------------------------------------------------------
#  # Convert BNG to WSGS84
#  osg_parse(grid_refs = "SN853872", coord_system = "WGS84")

## -----------------------------------------------------------------------------
#  osg_parse(grid_refs = someStations$`grid-reference`$ngr)

## ---- fig.width=7-------------------------------------------------------------
#  # Fetch only time series data from the waterml2 service
#  info <- cmr(id = "3001")
#  plot(info)
#  
#  # Fetch time series data and metadata from the waterml2 service
#  info <- cmr(id = "3001", metadata = TRUE)
#  plot(info$data, main=paste("Monthly rainfall data for the",
#                             info$meta$stationName,"catchment"),
#       xlab="", ylab=info$meta$units)

## ---- fig.width=7-------------------------------------------------------------
#  # Fetch only time series data
#  info <- gdf(id = "3001")
#  plot(info)
#  
#  # Fetch time series data and metadata from the waterml2 service
#  info <- gdf(id = "3001", metadata = TRUE)
#  plot(info$data, main=paste0("Daily flow data for the ",
#                              info$meta$station.name,
#                              " catchment (",
#                              info$meta$data.type.units, ")"))

## ---- fig.width=7-------------------------------------------------------------
#  # Search data/metadata
#  s <- cmr(c(3002,3003), metadata = TRUE)
#  
#  # s is a list of 2 objects (one object for each site)
#  plot(s[[1]]$data,
#       main = paste(s[[1]]$meta$station.name, "and", s[[2]]$meta$station.name))
#  lines(s[[2]]$data, col="green")

## -----------------------------------------------------------------------------
#  library(DT)
#  datatable(catalogue())

## -----------------------------------------------------------------------------
#  library(leaflet)
#  
#  leaflet(data = someStations) %>% addTiles() %>%
#    addMarkers(~longitude, ~latitude, popup = ~as.character(paste(id,name)))

## -----------------------------------------------------------------------------
#  library(dygraphs)
#  dygraph(info$data) %>% dyRangeSelector()

## -----------------------------------------------------------------------------
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
#  # Calculate the mean flow for each catchment
#  someStations$meangdf <- unlist(lapply(s2, mean))
#  
#  # Linear model
#  library(ggplot2)
#  ggplot(someStations, aes(x = as.numeric(`catchment-area`), y = meangdf)) +
#    geom_point() +
#    stat_smooth(method = "lm", col = "red") +
#    xlab(expression(paste("Catchment area [Km^2]",sep=""))) +
#    ylab(expression(paste("Mean flow [m^3/s]",sep="")))

