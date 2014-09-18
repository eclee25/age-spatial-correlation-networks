
## Name: Elizabeth Lee
## Date: 9/15/14
## Function: Create "nodes.csv" for jflowmap implementation from existing lat/long file. Zip3 (code), Zip3 (name), Lat, Lon
## Lat/Lon coordinates need to be converted to x,y coordinates in the correct map projection (designated by the map image/shapefile used in jflowmaps)

## Filenames: /home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/explore/mapping_code/Coord3digits.csv; /home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/explore/R_export/allpopstat_zip3_season_cl.csv (9/15/14: older version)
## Data Source: 
## Notes: https://code.google.com/p/jflowmap/wiki/HowToPrepareData
# us-background.png (from Eamon) uses the 'azequalarea' map projection; converted coordinates should be multipled by 100

## 
## useful commands:
## install.packages("pkg", dependencies=TRUE, lib="/usr/local/lib/R/site-library") # in sudo R
## update.packages(lib.loc = "/usr/local/lib/R/site-library")


dfsumm<-function(x) {
	if(!class(x)[1]%in%c("data.frame","matrix"))
		stop("You can't use dfsumm on ",class(x)," objects!")
	cat("\n",nrow(x),"rows and",ncol(x),"columns")
	cat("\n",nrow(unique(x)),"unique rows\n")
	s<-matrix(NA,nrow=6,ncol=ncol(x))
	for(i in 1:ncol(x)) {
		iclass<-class(x[,i])[1]
		s[1,i]<-paste(class(x[,i]),collapse=" ")
		y<-x[,i]
		yc<-na.omit(y)
		if(iclass%in%c("factor","ordered"))
			s[2:3,i]<-levels(yc)[c(1,length(levels(yc)))] else
		if(iclass=="numeric")
			s[2:3,i]<-as.character(signif(c(min(yc),max(yc)),3)) else
		if(iclass=="logical")
			s[2:3,i]<-as.logical(c(min(yc),max(yc))) else
			s[2:3,i]<-as.character(c(min(yc),max(yc)))
			s[4,i]<-length(unique(yc))
			s[5,i]<-sum(is.na(y))
			s[6,i]<-!is.unsorted(yc)
	}
	s<-as.data.frame(s)
	rownames(s)<-c("Class","Minimum","Maximum","Unique (excld. NA)","Missing values","Sorted")
	colnames(s)<-colnames(x)
	print(s)
} 

# outside function
# http://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

require(mapproj)

# read lat/long data
setwd('/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/outside_source_data')
orig <- read.csv('Coord3digits.csv', colClasses='character', header=T, na.strings='')
# read unique zip3 
setwd('/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/SQL_export')
uq_zips <- read.csv('unique_zip3.csv', colClasses='character', header=T)

# remove Coord3 data from outside of continental US
excludeSTATE <- c('PR', 'VI', 'HI', 'AK', 'AS', 'FM', 'GU', 'MH', 'MP', 'PW', 'AA', 'AP', 'AA')
orig2 <- orig[!(orig$STATE %in% excludeSTATE),]
new <- data.frame(zip3=orig2$zip3, lat=as.numeric(orig2$lat), long=as.numeric(orig2$long), stringsAsFactors=FALSE)
# replace na values with 0 in data frame
new[is.na(new)] <- 0
# add leading zeros to zip3
new$zip3 <- paste("00", new$zip3, sep='')
new$zip3 <- substrRight(new$zip3, 3)
# remove duplicates of zip3
new3 <- new[!duplicated(new$zip3),]
# include only zip3s in SDI dataset
new4 <- new3[(new3$zip %in% uq_zips$PATIENT_ZIP3),]
# include only zip3s where lat/long is not 0,0
new4 <- new4[!(new4$lat==0),]
new4 <- new4[!(new4$long==0),] # none were removed in this step

# perform checks that only zip3s in SDI dataset are in final dataset
sum(new4$zip %in% uq_zips$PATIENT_ZIP3) # 903
length(new4$zip %in% uq_zips$PATIENT_ZIP3) # 903
sum(uq_zips$PATIENT_ZIP3 %in% new4$zip) # 903
length(uq_zips$PATIENT_ZIP3 %in% new4$zip) #903

# jflowmap formatting #
# convert latlongs to correct projection coordinates 
# us-background.png (from Eamon) uses the 'azequalarea' map projection; converted coordinates should be multipled by 100
projection <- "azequalarea"
orientation <- c(30.82,-98.57,0)
# scale of 1 doesn't allow for good edge-bundling with this projection
scale <- 100
proj <- mapproject(y=as.numeric(new4$lat), x=as.numeric(new4$long), projection=projection, orientation=orientation)

# export flows projected coordinates
exp <- data.frame(Code=as.character(new4$zip3), Name=as.character(new4$zip3), stringsAsFactors=FALSE)
exp$Lat <- proj$y*scale
exp$Lon <- proj$x*scale

# 9/17/14 13:25
setwd('/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/jflowmap/input_files')
write.csv(exp, 'zip_nodes.csv', row.names=FALSE)

# export lat long coordinates for gephi
exp2 <- data.frame(zip3=as.character(new4$zip3), lat=as.numeric(new4$lat), long=as.numeric(new4$long), stringsAsFactors=FALSE)

# 9/18/14 11:33
setwd('/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/gephi/input_files')
write.csv(exp2, 'zip3_latlong.csv', row.names=FALSE)

# export list of zip3s in territories, AK, and HI
terr_zips <- orig[(orig$STATE %in% excludeSTATE),]$zip3
terrz <- paste("00", terr_zips, sep='')
terrz2 <- substrRight(terrz, 3)
# 9/16/14 18:28
setwd('/home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/R_export')
write.csv(terrz2, 'non_continental_zip3s.csv', row.names=FALSE)
