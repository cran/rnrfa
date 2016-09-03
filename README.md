


<!-- Edit the README.Rmd only!!! The README.md is generated automatically from README.Rmd. -->

rnrfa: An R package to Retrieve, Filter and Visualize Data from the UK National River Flow Archive
---------------

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.61439.svg)](http://dx.doi.org/10.5281/zenodo.61439)
[![Build Status](https://travis-ci.org/cvitolo/rnrfa.svg)](https://travis-ci.org/cvitolo/rnrfa.svg?branch=master)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/cvitolo/rnrfa?branch=master&svg=true)](https://ci.appveyor.com/project/cvitolo/rnrfa)
[![codecov.io](https://codecov.io/github/cvitolo/rnrfa/coverage.svg?branch=master)](https://codecov.io/github/cvitolo/rnrfa?branch=master)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/rnrfa)](https://cran.r-project.org/package=rnrfa)
[![CRAN Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/rnrfa)](https://cran.r-project.org/package=rnrfa)
[![CRAN Monthly Downloads](http://cranlogs.r-pkg.org/badges/rnrfa)](https://cran.r-project.org/package=rnrfa)

The UK National River Flow Archive serves daily streamflow data, spatial rainfall averages and information regarding elevation, geology, land cover and FEH related catchment descriptors.

There is currently an API under development that in future should provide access to the following services: metadata catalogue, catalogue filters based on a geographical bounding-box, catalogue filters based on metadata entries, gauged daily data for about 400 stations available in WaterML2 format, the OGC standard used to describe hydrological time series.  

The information returned by the first three services is in JSON format, while the last one is an XML variant.

The RNRFA package aims to achieve a simpler and more efficient access to data by providing wrapper functions to send HTTP requests and interpret XML/JSON responses.

### Dependencies
The rnrfa package depends on a number of CRAN packages. Check for missing dependencies and install them:

The rnrfa package is dependent on the **gdal** library and a number of CRAN packages. In unix-based operating systems **gdal** can be installed running the following command in terminal: 

    sudo apt-get install -y r-cran-rgdal


**R package dependencies** can be installed running the following code:


```r
install.packages(c("cowplot", "plyr", "httr", "xml2", "stringr", "xts", "rjson", "ggmap", "ggplot2", "sp", "rgdal", "parallel"))
```

This demo makes also use of external libraries. To install and load them run the following commands:


```r
packs <- c("devtools", "DT", "leaflet")
install.packages(packs)
lapply(packs, require, character.only = TRUE)
```

### Installation

The stable version of the **rnrfa** package is available from CRAN:


```r
install.packages("rnrfa")
```

Or you can install the development version from Github with [devtools](https://github.com/hadley/devtools):


```r
devtools::install_github("cvitolo/rnrfa")
```

Now, load the rnrfa package:


```r
library(rnrfa)
```

### Functions

#### List of monitoring stations
The R function that deals with the NRFA catalogue to retrieve the full list of monitoring stations is called catalogue(). The function, used with no inputs, requests the full list of gauging stations with associated metadata. The output is a dataframe containing one record for each station and as many columns as the number of metadata entries available. 


```r
# Retrieve information for all the stations in the catalogue:
allStations <- catalogue()
```

Those entries are briefly described as follows:

* `id` = Station identification number
* `name` = Name of the station
* `location` = Area in which the station is located
* `river` = River catchment
* `stationDescription` = General station description, containing information on weirs, ratings, etc.
* `catchmentDescription` = Information on topography, geology, land cover, etc.
* `hydrometricArea` = UK hydrometric area identification number
* `operator` = UK measuring authorities
* `haName` = Hydrometric Area name
* `gridReference` = OS Grid Reference number
* `stationType` = Type of station (e.g. flume, weir, etc.)
* `catchmentArea` = Catchment area in (Km^2)
* `gdfStart` = Year in which recordings started
* `gdfEnd` = Year in which recordings ended
* `farText` = Information on the regime (e.g. natural, regulated, etc.)
* `categories` = various tags (e.g. FEH\_POOLING, FEH\_QMED, HIFLOWS\_INCLUDED)
* `altitude` = Altitude measured in metres above Ordnance Datum or, in Northern Ireland, Malin Head.
* `sensitivity` = Sensitivity index calculated as the percentage change in flow associated with a 10 mm increase in stage at the $Q_{95}$ flow.
* `lat` = a numeric vector of latitude coordinates.
* `lon` = a numeric vector of longitude coordinates.

#### Station filtering
The same function catalogue() can be used to filter stations based on a bounding box or any of the metadata entries. 


```r
# Define a bounding box:
bbox <- list(lonMin=-3.82, lonMax=-3.63, latMin=52.43, latMax=52.52)

# Filter stations based on bounding box
someStations <- catalogue(bbox)
                                  
# Filter stations belonging to a certain hydrometric area
someStations <- catalogue(columnName="haName", columnValue="Wye (Hereford)")

# Filter based on bounding box & metadata strings
someStations <- catalogue(bbox,
                          columnName="haName",
                          columnValue="Wye (Hereford)")

# Filter stations based on threshold
someStations <- catalogue(bbox,
                          columnName="catchmentArea",
                          columnValue=">1")

# Filter based on minimum recording years
someStations <- catalogue(bbox,
                          columnName="catchmentArea",
                          columnValue=">1",
                          minRec=30)
                                  
# Filter stations based on identification number
someStations <- catalogue(columnName="id",
                          columnValue=c(3001,3002,3003))
                               
# Other combined filtering
someStations <- catalogue(bbox,
                          columnName="id",
                          columnValue=c(54022,54090,54091,54092,54097),
                          minRec=35)
```

#### Conversions
The only geospatial information contained in the list of station in the catalogue is the OS grid reference (column "gridRef"). The RNRFA package allows convenient conversion to more standard coordinate systems. The function "osg_parse()", for example, converts the string to easting and northing in the BNG coordinate system (EPSG code: 27700), as in the example below:


```r
# Where is the first catchment located?
someStations$gridReference[1]

# Convert OS Grid reference to BNG
osg_parse("SN853872")
```

The same function can also convert from BNG to latitude and longitude in the WSGS84 coordinate system (EPSG code: 4326) as in the example below.


```r
# Convert BNG to WSGS84
osg_parse("SN853872", CoordSystem = "WGS84")
```

osg_parse() also works with multiple references:


```r
osg_parse(someStations$gridReference)
```

#### Get time series data

The first column of the table "someStations" contains the id number. This can be used to retrieve time series data and convert waterml2 files to time series object (of class zoo). 

The National River Flow Archive serves two types of time series data: gauged daily flow and catchment mean rainfall.

These time series can be obtained using the functions gdf() and cmr(), respectively. Both functions accept three inputs: 

  * `id`, the station identification numbers (single string or character vector).

  * `metadata`, a logical variable (FALSE by default). If metadata is TRUE means that the result for a single station is a list with two elements: data (the time series) and meta (metadata).

  * `cl`, This is a cluster object, created by the parallel package. This is set to NULL by default, which sends sequential calls to the server.

Here is how to retrieve mean rainfall (monthly) data for _Shin at Lairg (id = 3001)_ catchment.


```r
# Fetch only time series data from the waterml2 service
info <- cmr(id = "3001")
plot(info)

# Fetch time series data and metadata from the waterml2 service
info <- cmr(id = "3001", metadata = TRUE)
plot(info$data, main=paste("Monthly rainfall data for the",
                           info$meta$stationName,"catchment"), 
     xlab="", ylab=info$meta$units)
```

Here is how to retrieve (daily) flow data for _Shin at Lairg (id = 3001)_ catchment.


```r
# Fetch only time series data from the waterml2 service
info <- gdf(id = "3001")
plot(info)

# Fetch time series data and metadata from the waterml2 service
info <- gdf(id = "3001", metadata = TRUE)
plot(info$data, main=paste("Daily flow data for the",
                           info$meta$stationName,"catchment"), 
     xlab="", ylab=info$meta$units)
```

#### Multiple sites
By default, the functions `getTS()` can be used to fetch time series data from multiple site in a sequential mode (using 1 core):


```r
# Search data/metadata in the waterml2 service
s <- cmr(c(3002,3003), metadata = TRUE)

# s is a list of 2 objects (one object for each site)
plot(s[[1]]$data, 
     main = paste(s[[1]]$meta$stationName, "and", s[[2]]$meta$stationName))
lines(s[[2]]$data, col="green")
```

### Interoperability

Upgrade your data.frame to a data.table:


```r
library(DT)
datatable(catalogue(all=FALSE))
```

Create interactive maps using leaflet:


```r
library(leaflet)

leaflet(data = someStations) %>% addTiles() %>%
  addMarkers(~lon, ~lat, popup = ~as.character(paste(id,name)))
```

Interactive plots using dygraphs:


```r
library(dygraphs)
dygraph(info$data) %>% dyRangeSelector()
```

Sequential vs Concurrent requests: a simple benchmark test

```r
library(parallel)
# Use detectCores() to find out many cores are available on your machine
cl <- makeCluster(getOption("cl.cores", detectCores()))

# Filter all the stations within the above bounding box
someStations <- catalogue(bbox)

# Get flow data with a sequential approach
system.time( s1 <- gdf(someStations$id, cl = NULL) )

# Get flow data with a concurrent approach (using `parLapply()`)
system.time( s2 <- gdf(id = someStations$id, cl = cl) )  
```

The measured flows are expected to increase with the catchment area. Let's show this simple regression on a plot:


```r
# Calculate the mean flow for each catchment
someStations$meangdf <- unlist( lapply(s2, mean) )

# Linear model
library(ggplot2)
ggplot(someStations, aes(x = as.numeric(catchmentArea), y = meangdf)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  xlab(expression(paste("Catchment area [Km^2]",sep=""))) + 
  ylab(expression(paste("Mean flow [m^3/s]",sep="")))
```

### Terms and Conditions
Please refer to the following Terms and Conditions for use of NRFA Data and disclaimer: http://nrfa.ceh.ac.uk/costs-terms-and-conditions

This package uses a non-public API which is likely to change. Package and functions herein are provided as is, without any guarantee.

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

### Meta

* Please [report any issues or bugs](https://github.com/cvitolo/rnrfa/issues).
* License: [GPL-3](https://opensource.org/licenses/GPL-3.0)
* Get citation information for `rnrfa` in R doing `citation(package = 'rnrfa')`
