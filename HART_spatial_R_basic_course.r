###code for Chris Muller class, a basic intro to GIS


#TODO:
  # fix my rgdal problem!
  # attend to all NOTES and TODOS in the body, and leave only those intended for students
  # clean up and finalize
  # add one ggplot2 map at the end of section 1.)




######################
# SECTION 0.) PREPARATION


  #Some of you may be brand spankin' new to R
  #If that's the case, today might be a whirlwind for you, but please ask questions! (No shame.)
  #I will do my best to explain things as I go, while also keeping time constraints in mind
  #Also, here's a good starting tutorial from our Stats Dept: http://www.stat.berkeley.edu/~spector/R.pdf


#Some good, boilerplate code for getting everyone's workspace, working directory, and packages set up
#(Gleaned from Shinhye Choi's and Patty Frontiera's spatial-R workshop, which you should register for if you'd
#be interested in learning more after today's material! Check http://gif.berkeley.edu/support/workshops.html
#for scheduling of this and other workshops.)

#clear workspace
rm(list = ls())

#set your working directory
setwd("/home/ihavehands/Hidden_Desktop/berk/dlab/consulting/Chris_Muller/") #NOTE: REPLACE WITH YOUR LOCAL DIRECTORY

#check for packages, install those you don't already have
required.pkg <- c('utils', 'raster', 'rgdal', 'sp', 'maptools', 'readxl', 'dismo')   #NOTE: LIST ALL PACKAGES HERE
pkgs.not.installed <- required.pkg[!sapply(required.pkg, function(p) require(p, character.only=T))]
install.packages(pkgs.not.installed, dependencies=TRUE)

# Now load all libraries 
lapply(required.pkg, library, character.only = TRUE)         





#Download the data for today's examples
  #Be sure to download into whatever directory you'd like to use as your current working directory
  #In case the following lines give you trouble, you can always just click this URL and then unzip the downloaded file manually:
    # https://www.dropbox.com/sh/avo94zu7nijgxw2/AAAkFrVMIeOTwGjnrhwsXjHVa?dl=1
download.file('https://www.dropbox.com/sh/avo94zu7nijgxw2/AAAkFrVMIeOTwGjnrhwsXjHVa?dl=0', 'example_files.zip', method = 'auto')
unzip('example_files.zip')









####################
# SECTION 1.) MAPPING DELL (2010)

#NOTE: PLAN THIS OUT

#GEOSPATIAL DATA is any data, whatever the form, that includes data associating it with its proper location on
#the face of the earth.

#Typically, this is expressed using COORDINATE PAIRS indicating the position of each point in the dataset. 

#Coordinates can be:
  #GEOGRAPHIC --> expressed as angles in two dimensions that determine points on a globe, usually using LATITUDE and LONGITUDE, 
  #CARTOGRAPHIC --> expressed as points on Euclidean (planar) axes, after the globe's surface has been PROJECTED onto the plane

#A PROJECTION is a scheme for representing the globe's surface on a flat plane
  #You can literally imagine placing the global inside a cylinder (a cylindrical projection) or a cone
  #(conical projection), putting a light inside it (jack-o-lantern-style), capturing the image on the cylinder
  #or cone, and then cutting it and unfolfing it flat
  #There are a number of other schemes
  #Importantly, it is proven impossible to devise a projection that maintains the true areas, sizes, angular
  #relations, and shapes of geographic entities (something's gotta give), so when choosing projections we often
  #want to consider which of these are most important to preserve for our purposes (i.e. do we want an 
  #equal-area projection? an equal-angle projection? etc.)


#Geospatial data typically is stored in one of two common ways:

  #VECTOR --> points, lines, or polygons (or sets of these) can be expressed as connected series of points;
  #the number of points used per real-world distance will determine the resolution of this representation (and
  #thus, the level or grain of detail with which we are representing the real-world complexity of a geographic
  #feature (consider, for example, that the most accurate map of the CA coast would be a 1:1 scaled replica of
  #the coast...). Basically, think of a 'connect-the-dots' model.

  #RASTER --> continuously spatially distributed variables are often represented by gridded cells, each having 
  #a location (which can be canonically as the coordinates of its center, or of its lower-left corner, etcetera) 
  #and at least one value (for the variable of interest). Cells will have a fixed cell-size (typically
  #expressed as the real-world distance represented by a cell-side, in either distance or degrees), and this
  #cell-size will determine the resolution (i.e. level of detail) with which we represent the real-world
  #complexity of this variable. Basically, think of a square-grid 'color-by-number' model. (Indeed, digital
  #photos are saved as rasters, usually with each cell having 3 values, for red, green, and blue.)



#Okay, that was a VERY brief run-down. This should all make more sense as we start to play with actual data.
#Most of the data we will be working with is vector data, but we'll also briefly see some raster data at the
#end.  Here goes!...





#First, let's read in the shapefiles we want (modern Peru border, and Mita region border)

  #NOTE: WE CAN EITHER USE THE RGDAL PACKAGE (WHICH IS EASIER AND PREFERABLE, BUT RIGHT 
  #NOW I AM HAVING AN OPAQUE ERROR THAT I CAN'T SOLVE WHEN I TRY TO LOAD RGDAL), 
  #OR WE CAN USE THE MAPTOOLS PACKAGE, BUT THIS PACKAGE CAN'T READ THE SHAPEFILE'S .PRJ FILE,
  #SO WE WILL HAVE TO STIPULATE AND FEED IT THE PROJECTION/CRS AS A 'PROJ4STRING'
  #GO TO HTTP://SPATIALREFERENCE.ORG FOR HELP CONSTRUCTING THESE STRINGS


  #read in the Mita border (they are lines)
  #proj4 = "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  #proj4 = CRS("+proj=longlat + ellps=WGS84")
  mita_and_peru_proj4 = CRS("+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

  mita = readShapeLines('./Dell_raw_data_files/MitaBoundary.shp', proj4string = mita_and_peru_proj4)
  
  
  #and read in the Peru map (it is a polygon set)
  peru = readShapePoly('./Dell_raw_data_files/peru_nw.shp', proj4string = mita_and_peru_proj4)
  
  




#What do these things look like? What's in them? Let's poke around a bit:
  mita

  mita$ # <Tab>
      #This returns the column names in the 'SpatialLinesDataFrame'; $ is a common symbol for getting columns
      #from a data.frame (more on the data.frame class later...)

  mita@ # <Tab>
      #This returns the components of an S4 object
        #R has both S3 and S4 objects (and R5 objects)
        #'objects' are instantiations of 'classes', which are 
        #canonically defined data-structures used to represent a certain 'real thing' (e.g. PDF, animal, shapefile, etc.)

  mita@lines #get the lines component

  typeof(mita@lines) #what is this sub-structure?

  mita@lines[[1]] #double-bracket for subsetting a list

  typeof(mita@lines[[1]]) #and what is this component?

  mita@lines[[1]]@ # <Tab> again

  #And so on and so forth... You can see how this is a highly structured data object, used to store 
  #all of the bits and pieces we need in order to characterize spatial lines

  #The method we've been using is an excellent way to 'dig into' the structure and organization of any
  #object that you're reading in with any R package, in order to start making sense of how it's arranged
  #and what lives where. Good to get in the habit.
  #You can of course also read the docs, but often they're long, complicated, and make your eyes glaze over
  #instantly. This is a bit more interactive and intuitive.
  
  
  #Now, try to get down to the actual coordinate pairs, and extract the first 100:

#~~~~~~~~~~~~~~~
#YOUR CODE HERE:




#~~~~~~~~~~~~~~~











  #Here's the answer
  mita@lines[[1]]@Lines[[1]]@coords[1:100,]

  #And as it turns out, there is a method built in the return the coordinates for an object, obviating the
  #complicated syntax that we just typed. Here's the answer using that instead:
  coordinates(mita)[[1]][[1]][1:100,]


  #Let's plot these:
  plot(mita@lines[[1]]@Lines[[1]]@coords[1:100,]) #as points
  #or 
  plot(mita@lines[[1]]@Lines[[1]]@coords[1:100,], type = 'l') # as lines



  #Also, some of the lines' coordinate arrays have the same coordinate in the first and last row. 
  #What does this mean?
  plot(peru@lines[[1]]@Lines[[819]]@coords, type = 'l')
  #Makes sense!



  #This should begin to give you a sense of how vector data is organized.
  #You could easily imagine how you could use this approach to save lines and polylines, 
  #polygons, and of course points








#Now let's map these shapefiles:

  #quick-and-dirty map (border in black, Mita region in red)
  plot(peru)
  plot(mita, col = 'red')

  #Wait! They didn't plot together. Annoying. Try instead:
  plot(peru)
  lines(mita, col = 'red')




  

#Now let's add another data-source: the point locations of all district capitals

#Read in XLS file (district capitals) as data.frame, using the readxl package
  caps = read_excel('./Dell_raw_data_files/locations.xlsx')

  #Remember that we saw a SpatialLinesDataFrame before? And mentioned the data.frame class?
  #A data.frame can be thought of in these ways:
    
    #The programmatic equivalent of spreadsheet (the obvious way; the intended use)
    
    #An array in which different datatypes are allowed in each column (the inutitive way, but not actually how it's implemented)

    #A list of data vectors, with each vector (i.e. column) having a name and being constrained to the same length
    
  #The latter is in fact how the data.frame class is implemented:
  typeof(caps)

  #However, because we used the readxl package to read it in, it was read in as a derived class defined by
  #that package (tbl.df)
  #Let's convert it to just a plan ol' data.frame (to avoid problems downstream)
  caps = as.data.frame(caps)

  #Data frames are sort of the bailiwick of R, and have a ton of methods (functions that inhere to an object
  #based on its class) and other functionalities (especially through other packages) that will make your data
  #work all-powerful. Definitely spend some time learning more about data frames.


  


#Plot our points on top of our lines and polygons  
  #We were able to read these points in as a data frame like this, without any fancy packages, because R is
  #currently unaware that these are geospatial data. In other words, our lat and lon columns are just two of a
  #number of columns of data, and they just happen to be numeric data. Clearly, thi wouldn't work (or would be
  #very messy and complicated) with line or polygon data, because in each row, instead of a value for lat and
  #a value for lon, we would need a vector (or a vector of vectors of vectors...) of points! 

  #For points, storage in XLS (or CSV, more simply) is typical. However, we can't just read in a run, as this
  #would create a problem:

  plot(peru)
  points(caps$lon, caps$lat)  #NOTE: Remember that lon is x, lat is y!

  plot(caps$lon, caps$lat)
  lines(peru)

  #Where the heck are our points?
  #They're way the heck out there...
  x_min = min(min(caps$lon), min(sapply(peru@lines[[1]]@Lines, function(x)min(x@coords[,1]))))-100
  y_min = min(min(caps$lat), min(sapply(peru@lines[[1]]@Lines, function(x)min(x@coords[,2]))))-100
  x_max = max(max(caps$lon), max(sapply(peru@lines[[1]]@Lines, function(x)max(x@coords[,1]))))+100
  y_max = max(max(caps$lat), max(sapply(peru@lines[[1]]@Lines, function(x)max(x@coords[,2]))))+100
  plot(peru, xlim = c(x_min, x_max), ylim = c(y_min, y_max))
  points(caps$lon, caps$lat, col = 'red')
  
  
  
  #This is because we need give our points a projection, and that projection needs to be matched to the 
  #projection of our lines and polygons!
  #One way (there are often various) that we can do this is by using the sp package's SpatialPointsDataFrame function
  #to create a spatial points object, and then assigning it the same CRS as above:

    #First, create a proj4string for the caps' coordinates
    #NOTE: they are simply unprojected, i.e. geographic, coordinates
    #We will use the most common standard, the World Geodetic System's most updated ellipsoid from 1984 (hence 'WGS84') 
    #This is a geodetic datum and coordinate reference system (CRS), which is essentially a mathematical model
    #of the globe as an ellipsoid and a set of reference points against which any point on the ellipsoid's
    #surface can be located 
    caps_proj4 = CRS("+proj=longlat + ellps=WGS84")
    sp_caps = SpatialPointsDataFrame(cbind(caps$lon, caps$lat), data = caps, proj4string = caps_proj4) 
    
    
    #What did this do? Take a look at all of the information output from the following commands:
    class(sp_caps)

    str(sp_caps)

    summary(sp_caps)


   
    #Now, we will project this using the projection used for the Mita lines and Peru polygons
    #A projection is a means of representing spatial locations from a non-Euclidean surface (i.e. curved
    #plane) on a Euclidean surface
    #You can literally think of this as projecting the surface of the globe onto a plane (and indeed, there is
    #a diversity of ways this can be done)
    #We can do this digitally using established algorithms that we can access through existing packages
    #The rgdal package, a mainstay of the R spatial world, will do this for us. 

    #NOTE: I CAN'T GET THIS TO WORK UNTIL I FIGURE OUT MY ANNOYING PROBLEM THAT WON'T LET ME LOAD RGDAL!
    proj_caps = spTransform(sp_caps, mita_and_peru_proj4)

  
    
 

    #Now our map should plot!
    plot(peru)
    lines(mita, col = 'red')
    points(sp_caps, col = 'blue')



    
    #lastly, let's plot this again, adding some nice bells and whistles, and write the output to a PDF

    #start the graphics device driver for producing PDF plots, which will save our map to the filename we've provided here
    pdf('Dell_mapping_output.pdf') 
    #TODO: FINISH THIS UP, AFTER DEBUGGING MY rgdal ISSUE; PERHAPS MAKE A NICE ggplot2 MAP, AND BRIEFLY DISCUSS ggplot?
    
    
   
    #once done plotting, turn off the graphics device
    dev.off() 


    #If you're interested in making maps like this for publication, I definitely recommend getting acquainted
    #with ggplot2. Picking this code apart and toying with it should be a good place to start! 













#####################
# SECTION 2.) AFRICAN ETHNIC GROUPS 


#Load in and explore our data
 
  #first we need to create the proper proj4string for our africa datasets                                                                             #it makes our lives easier that they all are already projected identically, so we just need one string to rule them all...                   
  africa_proj4 = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
   
  #now load in the explorer data
  routes = readShapeLines('./Pre_Colonial_Africa_Explorer_Routes/Explorer_Routes_Final.shp', proj4string = africa_proj4)

  #load in the Africa data
  africa = readShapePoly('./Pre_Colonial_Africa_Explorer_Routes/African_base.shp', proj4string = africa_proj4)

  #and lastly, load the tribal groups data
  tribes = readShapePoly('./Murdock_shapefile/borders_tribes.shp', proj4string = africa_proj4)




  #now we should have a good, solid idea of the structure and content of these objects, based on the output
  #of:

  routes

  africa





#Subset the routes for only those corresponding to voyages in or before 1875

  #Based on that understanding, and our knowledge of some basic R subsetting operations, this should also be straightforward

  #What column of our SpatialLinesDataFrame will we operate on?
  routes$ # <Tab>

  #Looks like it makes the most sense to use the Year_End column, so that we will include only voyages that
  #were completed by the end of calendar year 1875

  use_routes = routes[routes$Year_End <= 1875, ] #Because this is a spatial data.frame, we can use the same basic notation we would use for subsetting a data.frame, but we can still plot it with proper spatial referencing once doing that; in other words, the Spatial*DataFrame calss is written to provide intuitive data.frame functionalities for spatial data objects




  #Now let's find the subset of tribal groups regions that were visited by pre-1875 voyages
    #i.e. in a geometric/programmatic sense, which tribal polygons intersect with the pre-1875 routes?
    #Here's a handy function to find exactly this!
    intersect_tribes = intersect(tribes, use_routes)

    #It appears that 180 of the 843 tribal group areas intersected with our pre-1875 routes
    intersect_tribes

    #NOTE: Seemingly contradictorily, the intersect function belongs to the raster package. However, it seems to serve
    #our purposes better because it is also able to assess the intersection between two vector (vis. non-raster) 
    #datasets, and return the result still in a SpatialPolygonsDataFrame, which is precisely what we want
    #The other obvious option would be gIntersection function, in the rgeos package, but its particular
    #functionality would actually be a little less convenient. Again, many ways to do a thing. In navigating
    #among them, Google will be your best friend. You'll want to search for your needs using the proper
    #jargon, since it serves a universal standard for discussing this stuff. You'll likely find nicely framed
    #questions and multiple (vote-ranked) answers for what you need, mainly on StackExchange (e.g.
    #http://gis.stackexchange.com/questions/157194/difference-between-intersect-raster-and-gintersection-rgeos-in-polygons-r)





#Add a raster data layer

    #We don't have a lot of time to work with raster data here, so we'll just briefly read some in and plot
    #it, to demonstrate what it looks like and give you a feel for how it differs from vector data

    #A quick way to get some very basic global datasets is with the getData function in the raster package.
    #We can retrieve climatic data (min temp, max temp, annual precip), elevation, admin boundaries, and
    #some other stuff. 

    #Let's get precip, at 10' resolution from WorldClim:

    prec = getData('worldclim', var = 'prec', res = 10)

    #Now plot it and we're done
    plot(prec)

    #Psych! What is this? It appears this gives us average monthly precip. But we want average annual.
    #Let's just sum them! 
    #To figure out how to do this, let's first introspect our structure:
    prec

    #We seem to have summary information presented to us very similarly to how it was for our vector datasets
    #This is the result of the years of standardizing development that have been put into R's core spatial
    #packages (thank you, stranger-people!) 
    #However, this object is an instantiation of a particular class designed for raster data.

    #The RasterStack is pretty straightforward, can be thought of exactly as it sounds. If each raster is a single layer of
    #gridded, valued cells, then this is just a stack of those layers. 12 layers, namely (i.e. monthly).

    #So let's sum the monthly values to get a single layer of average annual precip:
    yr_prec = sum(prec)
    #Truly. That's it. 


    #Now, because we have global data, we'll want to hone in on Africa. In other words, we'll want to feed in
    #a new bounding box. Why don't we just extract this from our largest vector layer?
    box= bbox(africa)
    plot(yr_prec, ext = box)

    #(Also, here's a neat alternative for getting a bounding box:
    plot(yr_prec)
    draw_box = drawExtent()
    plot(yr_prec, ext = draw_box)

    #Neat, eh?)



    #NOTE: I CAN'T CURRENTLY RUN THIS, BECAUSE OF MY ISSUES WITH RGDAL, BUT FOR SOME REASON IT DOESN'T SEEM TO
    #HAVE MATTER FOR MAPPING PURPOSES. NONETHELESS, I'M LEAVING THIS, AND IT SHOULD WORK ONCE I'VE DEBUGGED MY
    #RGDAL ISSUE
    #Now lastly, we'll need to reproject our raster to match our other data's projection:
    reproj_yr_prec = projectRaster(yr_prec, africa_proj4)

    #And now we should be ready to plot everything together!


    


#Now finally, to save our map to a PDF
#This will basically follow the same steps as above.
#Put it all together and give it a try! We want to do the following:
  # - Plot the average annual precip raster as a basemap
  # - Add the border of Africa and its tribal group areas
  # - Add the polygons of those areas that intersect the pre-1875 routes, in transparent grey
  # - Add the routes of the pre-1875 routes, with dotted red lines
  # - Save the map as 'africa_map_practice.pdf' in our current working directory

#~~~~~~~~~~~~~~~
#YOUR CODE HERE:




#~~~~~~~~~~~~~~~














#Answer (DON'T PEEK!):


pdf('africa_map_practice.pdf')
plot(yr_prec, ext = box)
plot(africa, add = T)
plot(tribes, add = T)
plot(intersect_tribes, col = '#45454555', add = T) 
#NOTE: the color is expressed as a hexadecimal string (string of digits between 0 and F), with the last
#digit-couple indicating transparency (values less than FF are increasingly transparent)
lines(use_routes, col = 'red', lwd = 2)  #lwd just makes the lines thicker, so easier to see
dev.off()


