# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Severe Weather Lab
# 
# -----
# 
# This lab is going to be a little bit different than the others we have done so far in Synoptic. The purpose is going to be to introduce you to more concepts and tools in Python while playing with severe weather data. One of the first things you will probably notice is that this assignment is not being given to you on a handout. Instead, I'm writing this laboratory in what is known as an [IPython Notebook](http://ipython.org/ipython-doc/dev/interactive/htmlnotebook.html). In summary, the IPython Notebook is "a living document". You can interactively write Python code and have output displayed directly in the browser (figures included). Additionally, you can intersperse text (similar to what you can do in Microsoft Word) throughout the document using the [Markdown](http://daringfireball.net/projects/markdown/) language. Although it may seem a bit "weird" at first, I can assure you that the IPython Notebook is extremely powerful. In fact, I develop almost 100% of my code in the IPython Notebook now.
# 
# This lab is broken down into several distinct pieces. Each piece is designed to introduce you to various concepts in Python, while also having you understand more about how severe weather plays out on the synoptic scale.
# 
# ------
# 
# ## Excercise 1: Plotting Severe Weather Reports
# 
# The Storm Prediction Center generates daily text files (specifically, Comma Separated Value files, or CSV files for short) of National Weather Service Forecast Office Preliminary Storm Reports. These reports are called preliminary because they are typically based on reports to the NWS and not based on actual storm surveys. After a severe weather event the NWS will attempt to verify the storm reports and around two months after an event will submit a final report to NWS Headquarters for inclusion in the publication [Storm Data](http://www.ncdc.noaa.gov/oa/climate/sd/), which is considered the official storm report database. The daily severe storm reports (for recent events) can be found at the [SPC Storm Reports Page](http://www.spc.noaa.gov/climo/online/#reports). The problem with the preliminary reports page is that they contain information on events that could be continuation of other events, or even erroneous events. The SPC also produces CSV files of the reports in Storm Data which can be found at [SPC WCM Data Page](http://www.spc.noaa.gov/wcm/#data).
# 
# I've gone ahead and downloaded the tornado files and combined them into a single large file for ease of use. This file should be located in the same directory as this IPython Notebook file. It's named **1950-2011_onetor.csv**. 

# <codecell>

# Import the modules needed for this excercise

import numpy    # Our math library; also used to read in the data
import matplotlib.pyplot as plt    # Our plotting library
from mpl_toolkits.basemap import Basemap    # Our map library
import datetime    # Our Date/Time library

# Put plots in the document
%pylab inline

# <codecell>

# Now read in the data
file_name = './1950_2011_onetor.csv'
data = np.recfromtxt(file_name, unpack=True, dtype=None, names=True, delimiter=',')

# <markdowncell>

# The variable "data" now has the contents of our file in unique arrays. We can find the variables of the array by the following

# <codecell>

for name in data.dtype.names:
    print name,

# <markdowncell>

# You can read [this PDF](http://www.spc.noaa.gov/wcm/data/SPC_severe_database_description.pdf) to understand what each of the variables actually are.
# 
# Now in order to plot a specific tornado track we need to make sure we have the following information:
# 
# + Beginning Point (consisting of begining latitude and longitude)
# + Ending Point (consisting of ending latitude and longitude)
# + Tornado Rating (if we want to only plot specific ratings)
# + Date
# 
# So let's put each one of those fields into it's own variable.

# <codecell>

date = data['date']
time = data['time']
rating = data['f']
slat = data['slat']
slon = data['slon']
elon = data['elon']
elat = data['elat']

# <markdowncell>

# Now, let's take a look at the date and time variables. In particular, we want to see what each value in the array looks like so we know how to determine which points to plot (if we wanted to only plot the tornado tracks for a specific time period).

# <codecell>

for i in range(date.shape[0]):
    if 1 == 0 :
        print date[i], time[i]

# <markdowncell>

# Now that we know the format of the data (a list of strings), we could write a lot of code to match data so we could know which dates to plot. Fortunately, a lot of that is done for us in the "[datetime](http://docs.python.org/2/library/datetime.html)" module. This module allows Python to create date objects that can be compared the same way numbers are compared. So what we need to do now is put the date and time strings into a datetime format.
# 
# We will be useing the "strptime" method of the datetime object. Notice how the method name is essentially the phrase "strip time". This method takes the a string and a "format" and converts the string into a datetime object based on the specified format. You can see a list of the formatters at the bottom of the datetime module page (linked above).

# <codecell>

dts = []
print date.shape
for i in range(date.shape[0]):
    datetime_object = datetime.datetime.strptime(date[i] + ' ' + time[i], '%m/%d/%y %H:%M:%S')
    dts.append(datetime_object)

# <markdowncell>

# So the variable "dts" now is a list of datetime objects. These objects now allow you to conduct math operations on dates. In other words:

# <codecell>

print dts[0]
print dts[100]
print dts[1]
# Check to see if a date/time is earlier than another date:
print dts[0] < dts[100]

# Check to see if two dates are equal:
print dts[0] == dts[1]

# What's the difference between two dates:
print dts[100] - dts[0]

# <markdowncell>

# This is incredibly powerful in terms of filtering reports based on dates...and here's how we can do it!

# <codecell>

#begin_date = datetime.datetime(1991, 4, 26, 12, 0, 0)    # format here is "year, month, date, hour, minute, second"
#end_date = datetime.datetime(1991, 4, 27, 12, 0, 0)
begin_date = datetime.datetime(1999, 5, 3, 0, 0, 0)    # format here is "year, month, date, hour, minute, second"
end_date = datetime.datetime(1999, 5, 4, 0, 0, 0)

valid_slats = []
valid_slons = []
valid_elats = []
valid_elons = []
valid_ratings = []
valid_dates = []

# Now we loop through our original data to pull out the points that are valid
for i in range(len(dts)):
    if dts[i] >= begin_date and dts[i] <= end_date:
        if rating[i] >= 4: # Filter for only violent tornadoes
           valid_dates.append(dts[i])
           valid_slats.append(slat[i])
           valid_slons.append(slon[i])
           valid_elats.append(elat[i])
           valid_elons.append(elon[i])
           valid_ratings.append(rating[i])
    
    
    

# <markdowncell>

# Now, let's plot them on a map! For this we are going to use the [Basemap](http://matplotlib.org/basemap/) module. Now, click on the Basemap link and you'll be taken to the documentation page for this module. (It's **incredibly** useful.) This "landing page" has links to the different steps necessary to create our map. In particular, the first place we need to start is "Setting Up The Map". Please click on this link. The text at the top of the page explains why we have to use a mapping module instead of just simply plotting. At the bottom, there are links to various different [map projections](http://en.wikipedia.org/wiki/Map_projection). I encourage you to look through the various map projections and find one that you like. For this lab, I'm going to use the [Lambert Conformal Projection](http://matplotlib.org/basemap/users/lcc.html), but you can eventually use any you would like, but for now, I strongly recommend that you stick with the Lambert Conformal Projection (LCC) for now.
# 
# Click on the LCC projection link and you'll be taken to a page that shows how to set up the map and gives an illustration of the map. In particular, the line (in the sample code) that actually constructs the map is the line that begins with "m = Basemap(...)". There are two different ways to construct the map:
# 
# 1. Giving a center point and then a width/height in meters
# 2. Giving the lower-left corner and upper-right corder latitudes and longitudes.
# 
# The example online shows the first way. I've used the second way below.
# 

# <codecell>

# setup lambert conformal basemap.
# resolution is the resolution of the map. 'c' is crude outlines, 'l' is low, 'i' is intermediate. 
#     as you increase resolution, it takes longer to draw, but you get more detail
# llcrnrlon is the longitude of the lower-left corner
# llcrnrlat is the latitude of the lower-left corner
# urcrnrlon is the longitude of the upper-right corner
# urcrnrlat is the latitude of teh upper-right corner
# lat_1 is first standard parallel
# lon_0 is the longitude that is perfectly north-south (determines how rotated your map looks)
# lat_2 is second standard parallel (defaults to lat_1).
# rsphere=(6378137.00,6356752.3142) specifies shape of the earth (WGS4 ellipsoid in this case)
# area_thresh=1000 means don't plot coastline features les
#     than 1000 km^2 in area.
if 1 == 1: 
    m = Basemap(resolution='i', projection='lcc', llcrnrlon=-120,
            llcrnrlat=22.5, urcrnrlon=-65, urcrnrlat=48.5, lat_1=38.5, lon_0=-97.5,
            lat_2=38.5, rsphere=(6378137.00, 6356752.3142), area_thresh=10000)
else:
    m = Basemap(width=6000000,height=5000000,
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=45.,lat_2=55,lat_0=45,lon_0=-95.)    
    
# The variable "m" is now our map. We can use this map instance to draw various things such as:
m.drawcoastlines()    # Draws the coastlines
m.drawcountries()     # Draw the country borders
m.drawstates()        # Draw the state borders
m.fillcontinents(color='#DADADA', lake_color='aqua')    # Fill the continents with a color
m.drawmapboundary(fill_color='aqua')     # This will plot the oceans the given color

# <markdowncell>

# Now we need to plot our tornado data! We can do that by:

# <codecell>

# Note, we don't have to recreate our map because we did it above.
# However, after we created the picture above, the map drawing 
# goes away, so we have to recreate that part

m.drawcoastlines()    # Draws the coastlines
m.drawcountries()     # Draw the country borders
m.drawstates()        # Draw the state borders
m.fillcontinents(color='#DADADA', lake_color='aqua')    # Fill the continents with a color
m.drawmapboundary(fill_color='aqua')     # This will plot the oceans the given color

# Now we loop through the tornado data to plot
for i in range(len(valid_dates)):
    x1 = valid_slons[i]
    y1 = valid_slats[i]
    x2 = valid_elons[i]
    y2 = valid_elats[i]
    m.plot([x1, x2], [y1, y2], 'rv', linestyle='solid', linewidth=1)

# <markdowncell>

# You'll notice that nothing happened. Nothing was plotted on our map! This is because the map is constructed in "projection space" and our tornado points are contructed in latitude and longitude space. In order to plot them on the map, we need to convert the latitude and longitudes into projection space. Fortunately, our map instance (m) will handle that conversion for us! Here's how:

# <codecell>

# Format is: new_x, new_y = m(old_x, old_y)

new_slons, new_slats = m(valid_slons, valid_slats)
new_elons, new_elats = m(valid_elons, valid_elats)

# <markdowncell>

# Now we can try our plot again, but this time using the projection points:

# <codecell>

# Note, we don't have to recreate our map because we did it above.
# However, after we created the picture above, the map drawing 
# goes away, so we have to recreate that part

m.drawcoastlines()    # Draws the coastlines
m.drawcountries()     # Draw the country borders
m.drawstates()        # Draw the state borders
m.fillcontinents(color='#DADADA', lake_color='aqua')    # Fill the continents with a color
m.drawmapboundary(fill_color='aqua')     # This will plot the oceans the given color

# Now we loop through the tornado data to plot
for i in range(len(valid_dates)):
    x1 = new_slons[i]
    y1 = new_slats[i]
    x2 = new_elons[i]
    y2 = new_elats[i]
    m.plot([x1, x2], [y1, y2], 'rv', linestyle='solid', linewidth=1)

# <markdowncell>

# Ugh, what happpend? Why are there so many tracks heading off to Europe?! Well, let's take a look at our latitudes and longitudes...

# <codecell>

for i in range(len(valid_slats)):
    if 1 == 0:
        print valid_slons[i], valid_slats[i], valid_elons[i], valid_elats[i]

# <markdowncell>

# You should notice that some of our ending latitudes and longitudes are 0. We'll need to fix this to plot correctly. The easiest solution is to set the ending points to the beginning points if the ending points are 0. So let's do that:

# <codecell>

for i in range(len(valid_slats)):
    if valid_elons[i] == 0 or valid_elats[i] == 0:
        valid_elons[i] = valid_slons[i]
        valid_elats[i] = valid_slats[i]

# <markdowncell>

# Now we have to convert this back into projection space:

# <codecell>

# Format is: new_x, new_y = m(old_x, old_y)

new_slons, new_slats = m(valid_slons, valid_slats)
new_elons, new_elats = m(valid_elons, valid_elats)

# <markdowncell>

# Now we plot again:

# <codecell>

# Note, we don't have to recreate our map because we did it above.
# However, after we created the picture above, the map drawing 
# goes away, so we have to recreate that part

m.drawcoastlines()    # Draws the coastlines
m.drawcountries()     # Draw the country borders
m.drawstates()        # Draw the state borders
m.fillcontinents(color='#DADADA', lake_color='aqua')    # Fill the continents with a color
m.drawmapboundary(fill_color='aqua')     # This will plot the oceans the given color

# Now we loop through the tornado data to plot
for i in range(len(valid_dates)):
    x1 = new_slons[i]
    y1 = new_slats[i]
    x2 = new_elons[i]
    y2 = new_elats[i]
    m.plot([x1, x2], [y1, y2], 'rv', linestyle='solid', linewidth=1)

# <markdowncell>

# Now, you'll notice that there aren't any tracks! This is because we are zoomed too far out on this map to see any tracks. Now, on your own, play around with the map settings [m = Basemap(...) line] to configure a map that is zoomed in to see tornado tracks. Use the code block below. Note, you may need to combine things from multiple code blocks above.

# <codecell>

m = Basemap(width=400000,height=400000,
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=45.,lat_2=55,lat_0=36,lon_0=-97.)    
m.drawcoastlines()    # Draws the coastlines
m.drawcountries()     # Draw the country borders
m.drawstates()        # Draw the state borders
m.fillcontinents(color='#DADADA', lake_color='aqua')    # Fill the continents with a color
m.drawmapboundary(fill_color='aqua')     # This will plot the oceans the given color

new_slons, new_slats = m(valid_slons, valid_slats)
new_elons, new_elats = m(valid_elons, valid_elats)
for i in range(len(valid_dates)):
    x1 = new_slons[i]
    y1 = new_slats[i]
    x2 = new_elons[i]
    y2 = new_elats[i]
    m.plot([x1, x2], [y1, y2], 'rv', linestyle='solid', linewidth=1)


# <markdowncell>

# Now, play around with the code until you understand it. Change variables and see what happens. Change the time period to an event that interests you! (This is how you learn!).
# 
# Once you begin to feel comfortable, do one of the following (or both!) in the code block below (Yes. This is an assignment):
# 
# 1. On your zoomed map, plot only the violent tornadoes (rating >= 4).
# 2. Color code your tornado symbols based on their rating
# 
# (Note: if you do both of these, you'll need to add an additional code block.)

# <codecell>


# <markdowncell>

# -----
# ## Plotting Environment Data for Your Event
# 
# Now that we have plotted the tornadoes for a given tornado event, let's plot some of the environment variables that go along with the tornado event. For this we will use the North American Regional Reanalysis (NARR) datasets available at the Earth System Research Laboratory (ESRL)'s Physical Science Division (PSD) [NARR plotting page](http://www.esrl.noaa.gov/psd/cgi-bin/data/narr/plothour.pl). On this page you can choose from a list of variables and levels, supply a date (back to January 1979) and then create plots of said variable. Take a minute to play around with the website and get comfortable creating plots.
# 
# One of the nice things about this website is that after you create a plot, you can download a netCDF file of the data you plotted. To do this, take a look on the lower-left side of the webpage that contains a plot you generated. You should see a link called: "Get a copy of the netcdf data file used for the plot". If you click on this link, you will download the netCDF file of the data you plotted. We will make significant use of this feature in this part of the lab. **I strongly recommend you change the name of the downloaded file to something more descriptive of what the file contains so you'll be able to quickly identify what is in which file.**
# 
# What I would like you to do is download netCDF files of fields you deem important and create plots of them in code blocks below. In each of the plots, I would like to see the tornado tracks also plotted. This way, when you are done, you will have maps of the synoptic scale environment and can relate them to where tornadoes actually occrred.  Fields that I would like to see plotted:
# 
# + 500 mb Heights
# + 500 hPa Winds
# + Pressure Reduced to Mean Sea Level (prmsl)
# + 10-m Winds
# + CAPE
# + 2-m Dewpoint
# + At least one field of your own choosing
# 
# Generally speaking, I am not going to specifically tell you what fields to put on which map. You can combine them however you see fit.
# 
# **After you finish creating your maps, please write up a 1-page synopsis of the event you chose.** This write up should basically explain the maps you've created using what you've learned from synoptic class. For examples of what I'm looking for, check out the [SPC Severe Weather Outlooks](http://www.spc.noaa.gov/products/outlook/day1otlk.html). Please do not just copy (or rewrite) the SPC discussion for your event. Trust me, I will check your write up and compare it to the SPC write up for the given event.
# 
# When you are finished, please email your IPython Notebook and your write up to **pmarshwx+synoptic@gmail.com**
# 
# -----
# ### Some Useful Tips 
# 
# + The matplotlib module has an excellent webpage. You can find all the plotting commands [here](http://matplotlib.org/api/pyplot_summary.html)
# + Also, make use of the [matplotlib gallery](http://matplotlib.org/gallery.html) to see sample plots along with the code to generate them.
# + To read in a netcdf file you will need to import the netCDF4 module. This can be done by: **import netCDF4 as nc**.
# + The commands to read in netcdf files is: **f = nc.Dataset(file_path_and_name, 'r')**.
# + To see what variables are contained in the file you can: **print f.variables**
# + To read in a variable from the netcdf file you: **variable = f.variables['variable name'][:]**
# > You must be careful here! Each variable can be of different dimensions (different levels, etc).
# > > + You'll need to read in the latitudes and longitudes from the file (don't forget to change them to projection space for plotting!). These will be two dimensional arrays. They can be read in by: **var = f.variable['var name'][:]**
# > > + The "field" variables will be four dimensional arrays, but you only need two of them. To read these variables in use: **var = f.variables['var name'][0,0,:,:]**
# 
# + To plot contours you will use the [contour](http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.contour) argument. Example: **plt.contour(x, y, var, levels=range(min, max+1, interval), colors='color name', linewidth=number, linestyle='solid')** See the contour link for more information.
# + To plot a color-filled contour you'll use the [contourf](http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.contourf) argument. It's used similarly as the contour argument.
# + **Note, when you are plotting on a map, replace the "plt" above with "m" (or whatever you called your map instance).**
# + To add a colorbar, use **plt.colorbar()**. If you have multiple color-filled instances on the same map, you will need to do the following: 
# 
# >> fill1 = plt.contourf(...)
# 
# >> fill2 = plt.contourf(...)
# 
# >> plt.colorbar(fill1)
# 
# >> plt.colorbar(fill2)
# 
# + Don't forget to add your tornado data! You can get that code from the previous section!
# 
# 
# Additional information will be provided should repeated questions arise. These "tips" will be posted in the **Club 4424 Facebook page**. If you are not a member of the group, please let me know.

# <codecell>


