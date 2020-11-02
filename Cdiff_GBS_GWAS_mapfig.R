#Cdiff GBS GWAS map figure
#8/18/17

library(rgdal) # Commands for reprojecting the vector data.
library(rworldmap) # Recently updated mapping program.
library(rworldxtra) # Add-ons for rworldmap.

#quick redo, color countries, two legend, no germ trial

# pdf("KTurnerFig1.pdf", useDingbats=FALSE, width=13.38)
# png("Cdiff_GWAS_mapfig.png", width=1200, height = 600, pointsize = 16)
#svg("collectionMap_bw.svg", pointsize = 12)
# setEPS( horizontal = FALSE, onefile = FALSE, paper = "special") #340 mm is c. 13 in, height and width in in
# postscript("colMap_bw.eps")

jpeg("fig1_rangeMap.jpg", width=1200, height = 600, pointsize = 10, quality = 100, res = 300)

# setEPS()
# postscript("figCdiff_rangeMap.eps", onefile = FALSE, paper = "special", height = 4.5, width = 7) #horizontal = FALSE, height = 7,

# eps <- function(file, onefile=FALSE, width=5.5, height = 5.5, horizontal=FALSE,paper="special", ...){
#   postscript(file=file,width=width,height=height,onefile=onefile,horizontal=horizontal,paper=paper,title=file,...)
#   par(bty='l')
# }

projectionCRS <- CRS("+proj=laea +lon_0=0.001 +lat_0=89.999 +ellps=sphere") #the ellps 'sphere' has a radius of 6370997.0m
par(mai=c(0,0,0.2,0)) #,xaxs="i",yaxs="i"
sPDF <- getMap()[-which(getMap()$ADMIN=='Antarctica')] 
sPDF <- spTransform(sPDF, CRS=projectionCRS)
setLims <- TRUE #FALSE back to whole world
#setLims <- FALSE
if ( !setLims )
{
  xlim <- ylim <- NA
} else
{
  ### TRY FIDDLING WITH THESE LIMITS ###
  xlimUnproj <- c(-52,120)
  ylimUnproj <- c(8,30)
  sPointsLims <- data.frame(x=xlimUnproj, y=ylimUnproj)
  coordinates(sPointsLims) = c("x", "y")
  proj4string(sPointsLims) <- CRS("+proj=longlat +ellps=WGS84")
  sPointsLims <- spTransform(sPointsLims, CRS=projectionCRS)
  xlim <- coordinates(sPointsLims)[,"x"]
  ylim <- coordinates(sPointsLims)[,"y"]  
}

# countries <- as.data.frame(cbind(sPDF$ADMIN))
# countries$country <- sPDF$ADMIN
# countries$color <- "lightgray"
# #countries <- reorder(countries$country, countries$V1)
# #countries <- unique(countries)
# 
# countries[countries$country %in% c("United States of America", "Canada"),]$color <- "red"
# grayco <- as.vector(countries$color)
# 
# mapCountryData(sPDF, nameColumnToPlot="SOVEREIGNT", mapTitle='Global range of C. diffusa', colourPalette=countries$color, borderCol ='gray24',addLegend = FALSE, xlim=xlim, ylim=ylim)

# sPDF <- getMap()
# #list of country names
# sPDF$ADMIN
#setup a color code column filled with 1's
sPDF$colCode <- 1
# tst <- sPDF@data
# View(tst)
#set codes for specified countries
sPDF$colCode[ which(sPDF$ADMIN %in% c("Canada","United States of America"))] <- 2
sPDF$colCode[ which(sPDF$ADMIN %in% c("Armenia","Azerbaijan", "Bulgaria", "Georgia", 
                                      "Greece", "Moldova", "Romania","Russia", "Turkey",
                                      "Ukraine", "Serbia"))] <- 3
sPDF$colCode[ which(sPDF$ADMIN %in% c("Poland", "Belarus", "Italy", "Syria", "Czech Republic",
                                      "Estonia", "Switzerland","Latvia","Lithuania", 
                                      "Slovenia", "Serbia","Austria","Belgium", "France",
                                      "Germany","Hungary","Luxembourg","Norway","Slovakia",
                                      "Spain", "United Kingdom", "Kazakhstan", "Turkmenistan", "China"))] <- 4
# tst <- sPDF@data
# View(tst)
#create a colour palette - note for each value not for each country
# colourPalette <- c("lightgray","#F8766D","#00BFC4", "cadetblue1")
# colourPalette <- c("lightgray","lightgray")
colourPalette <- c("lightgray","#666699","#66CC00", "#99FF99")
# spName <- plotmath(italic("Centaurea diffusa"))
par(mar=c(0,0,0,0))
mapCountryData(sPDF, nameColumnToPlot="colCode", mapTitle=NA,
               colourPalette=colourPalette, borderCol ='gray24', addLegend = FALSE,
               xlim=xlim, ylim=ylim, catMethod=c(0,1,2,3,4))
#note that catMethod defines the breaks and values go in a category if they are <= upper end
#mapTitle=bquote(Global~range~of~italic(Centaurea)~italic(diffusa)) 


# # to plot states you can get data from here :
# #   http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_1_states_provinces_shp.zip
# inHigh <- "ne_10m_admin_1_states_provinces_shp.shp"
# sPDFhigh <- readShapePoly(inHigh)
# # attributes <- sPDFhigh@data
# #test plotting
# mapPolys(sPDFhigh, nameColumnToPlot='admin', addLegend = FALSE)
# #then you could subset the countries you want, 
# sPDFhigh <- sPDFhigh[which(sPDFhigh$ADMIN %in% c("Canada","United States of America","Russia"))]
# # look up state names & do 
# mapPolys(sPDFhigh$name ,add=TRUE)
# # = c("Yukon","British Columbia", "Washington", "Oregon")

popcoord <- read.delim("CdiffGBS_GWAS_coord.txt", stringsAsFactors = F)
popcoord[5,4] <- "Assembled reference"
popcoord$pch <- 16 #invasive
popcoord[popcoord$Type=="Native C. diffusa",]$pch <- 17
popcoord[popcoord$Type=="Assembled reference",]$pch <- 15


coordinates(popcoord) = c("Longitude", "Latitude")
proj4string(popcoord) <- CRS("+proj=longlat +ellps=WGS84")
sPointsDF <- spTransform(popcoord, CRS=projectionCRS)
points(sPointsDF, pch=popcoord$pch, cex=0.75) #cex=1.5 for eps

llgridlines(sPDF, easts=c(-90,-180,0,90,180), norths=seq(0,90,by=15), 
            plotLabels=FALSE, ndiscr=1000) #ndiscr=num points in lines
#lat markings...

markings <- data.frame(Latitude=as.numeric(c(75,60,45,30,15,85,85)), Longitude=as.numeric(c(-45,-45,-45,-45,-45,0,180)),name=c("75", "60","45","30","15","0","180"))
coordinates(markings) = c("Longitude", "Latitude")
proj4string(markings) <- CRS("+proj=longlat +ellps=WGS84")
sPointsDFmark <- spTransform(markings, CRS=projectionCRS)
text(sPointsDFmark, labels = sPointsDFmark$name, cex=1) #pch2 for triangles #cex=1.2 for eps

# pole <- data.frame(x=0, y=90)
# coordinates(pole) = c("x", "y")
# proj4string(pole) <- CRS("+proj=longlat +ellps=WGS84")
# pole <- spTransform(pole, CRS=projectionCRS)
# points(pole, pch=8, cex=2, lwd=2)

pch_title <- expression(paste("Invasive ", italic("C. diffusa")))


legend("bottomleft", legend = c(expression(paste("Invasive ", italic("C. diffusa"))),expression(paste("Native ", italic("C. diffusa"))), "Assembled reference"),
       pch=c(16,17,15),  bg="white", title = "Sampled populations", cex=0.75) #cex=1 for eps
legend("topright", c("Invasive", "Native","Present, status unknown"), fill=c("#666699","#66CC00", "#99FF99"),
       title=expression(paste(italic("C. diffusa "), "range")), bg="white", cex=0.75) #cex=1 for eps
box(lty="solid", col = "black")
# #shameless plug !
# mtext("map made using rworldmap", line=-1, side=1, adj=1, cex=0.6)

# 
dev.off()