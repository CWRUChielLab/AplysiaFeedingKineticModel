import os
import sys
import csv
import matplotlib
#matplotlib.use('Agg') # force matplotlib to not use XWindows backend
from math import sqrt
from scipy.interpolate import spline
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.lines as mlines
from matplotlib.patches import Circle
from matplotlib import animation


######################################################
##                 Video Parameters                 ##
######################################################

#File Name
moviefile = "animation.mp4"
#size of movie
width = 1600
height = 400
#frames per second
fps = 100
#bitrate to reduce compression
br = 2000
#moviesize
dpi = 80


######################################################
##                    Read File                     ##
######################################################

#The following code reads the animationinfo.txt.txt file and converts it to a csv file
tabDelimitedFile = open("animationinfo.txt", "r")
splitFileToList = "animationinfo.txt".split(".")
listInput = []
for line in tabDelimitedFile:
    listLine = line.replace("\n","").split("\t")
    listInput.append(listLine)

tabDelimitedFile.close()
newFile = open(splitFileToList[0] + "-csvTranslated.csv", "w")
for line in listInput:
    writeLineForNewFile = ",".join(line)
    newFile.write(writeLineForNewFile + "\n")

newFile.close()

#The following code reads the animationinfo-csvTranslated.csv and saves the position in an array list "position"
position = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            position.insert(i,column[1])
            i+=1

#The following line prints the time array to confirm that it contains the correct values
#print(position)

#The following code reads the animationinfo-csvTranslated.csv and saves the radius in an array list "radius"
radius = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            radius.insert(i,column[2])
            i+=1

#The following line prints the time array to confirm that it contains the correct values
#print(radius)

#The following code reads the animationinfo-csvTranslated.csv and saves the angle in an array list "angle"
angle = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            angle.insert(i,column[3])
            i+=1

#The following line prints the time array to confirm that it contains the correct values
#print(angle)


#The following code reads the animationinfo-csvTranslated.csv and saves the xctop in an array list "xctop"
xctop = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            xctop.insert(i,column[4])
            i+=1

#The following line prints the xctop array to confirm that it contains the correct values
#print(xctop)


#The following code reads the animationinfo-csvTranslated.csv and saves the xcbottom in an array list "xcbottom"
xcbottom = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            xcbottom.insert(i,column[5])
            i+=1

#The following line prints the xcbottom array to confirm that it contains the correct values
#print(xcbottom)

#The following code reads the animationinfo-csvTranslated.csv and saves the ytop in an array list "ytop"
ytop = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            ytop.insert(i,column[6])
            i+=1

#The following line prints the ytop array to confirm that it contains the correct values
#print(ytop)

#The following code reads the animationinfo-csvTranslated.csv and saves the ybottom in an array list "ybottom"
ybottom = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            ybottom.insert(i,column[7])
            i+=1

#The following line prints the ybottom array to confirm that it contains the correct values
#print(ybottom)

#The following code reads the animationinfo-csvTranslated.csv and saves the y in an array list "i1i3radiusarray"
i1i3radiusarray = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            i1i3radiusarray.insert(i,column[8])
            i+=1

#The following line prints the i1i3radiusarray array to confirm that it contains the correct values
#print(i1i3radiusarray)

#The following code reads the animationinfo-csvTranslated.csv and saves the i2length in an array list "i2lengtharray"
i2lengtharray = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            i2lengtharray.insert(i,column[9])
            i+=1

#The following line prints the i2length array to confirm that it contains the correct values
#print(i2lengtharray)

#The following code reads the animationinfo-csvTranslated.csv and saves the topangle in an array list "topanglearray"
topanglearray = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            topanglearray.insert(i,column[10])
            i+=1

#The following line prints the topangle array to confirm that it contains the correct values
#print(topanglearray)

#The following code reads the animationinfo-csvTranslated.csv and saves the bottomangle in an array list "bottomanglearray"
bottomanglearray = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            bottomanglearray.insert(i,column[11])
            i+=1

#The following line prints the bottomangle array to confirm that it contains the correct values
#print(bottomanglearray)

#The following code reads the animationinfo-csvTranslated.csv and saves the furthestbackxpoint in an array list "furthestbackxpointarray"
furthestbackxpointarray = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            furthestbackxpointarray.insert(i,column[12])
            i+=1

#The following line prints the furthestbackxpointarray to confirm that it contains the correct values
#print(furthestbackxpointarray)

#The following code reads the animationinfo-csvTranslated.csv and saves the furthestbackypoint in an array list "furthestbackypointarray"
furthestbackypointarray = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            furthestbackypointarray.insert(i,column[13])
            i+=1

#The following line prints the furthestbackypointarray to confirm that it contains the correct values
#print(furthestbackypointarray)

#The following code reads the animationinfo-csvTranslated.csv and saves the i1i3contacttopy in an array list "i1i3contacttopyarray"
i1i3contacttopyarray = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            i1i3contacttopyarray.insert(i,column[14])
            i+=1

#The following line prints the i1i3contacttopyarray to confirm that it contains the correct values
#print(i1i3contacttopyarray)

#The following code reads the animationinfo-csvTranslated.csv and saves the i1i3contactbottomy in an array list "i1i3contactbottomyarray"
i1i3contactbottomyarray = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            i1i3contactbottomyarray.insert(i,column[15])
            i+=1

#The following line prints the i1i3contactbottomyarray to confirm that it contains the correct values
#print(i1i3contactbottomyarray)

#The following code reads the animationinfo-csvTranslated.csv and saves the ocontacttopx in an array list "ocontacttopxarray"
ocontacttopxarray = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            ocontacttopxarray.insert(i,column[16])
            i+=1

#The following line prints the ocontacttopxarray to confirm that it contains the correct values
#print(ocontacttopxarray)

#The following code reads the animationinfo-csvTranslated.csv and saves the ocontacttopy in an array list "ocontacttopyarray"
ocontacttopyarray = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            ocontacttopyarray.insert(i,column[17])
            i+=1

#The following line prints the ocontacttopyarray to confirm that it contains the correct values
#print(ocontacttopyarray)

#The following code reads the animationinfo-csvTranslated.csv and saves the ocontactbottomx in an array list "ocontactbottomxarray"
ocontactbottomxarray = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            ocontactbottomxarray.insert(i,column[18])
            i+=1

#The following line prints the ocontactbottomxarray to confirm that it contains the correct values
#print(ocontactbottomxarray)

#The following code reads the animationinfo-csvTranslated.csv and saves the ocontactbottomy in an array list "ocontactbottomyarray"
ocontactbottomyarray = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            ocontactbottomyarray.insert(i,column[19])
            i+=1

#The following line prints the ocontactbottomyarray to confirm that it contains the correct values
#print(ocontactbottomyarray)

#The following code reads the animationinfo-csvTranslated.csv and saves the bigxval in an array list "bigxvalarray"
bigxvalarray = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            bigxvalarray.insert(i,column[19])
            i+=1

#The following line prints the bigxvalarray to confirm that it contains the correct values
#print(bigxvalarray)

#The following code reads the animationinfo-csvTranslated.csv and saves the x1 in an array list "x1array"
x1array = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            x1array.insert(i,column[20])
            i+=1

#The following line prints the x1array to confirm that it contains the correct values
#print(x1array)


######################################################
##                Plot 2D Simulation                ##
######################################################

a = plt.subplot(111, aspect='equal')

#Define Shapes and Variables. discretemoment represents the first time step of the animationinfo-csvTranslated.csv file, and each time step becomes a frame in the final video.

def createShapes(discretemoment, figure):
    #Define Variables by acsessing array values
    
    oradius = float(radius[discretemoment]) #"a" in the kinetic model
    odiameter = 2*oradius #2*a
    ominoraxisradius = (5.0 * sqrt(5.0) / sqrt(oradius * 1000.0)) / 1000.0 #"b" in the kinetic model
    ominoraxisdiameter = 2*ominoraxisradius #2*b
    oangle = -(90 - float(angle[discretemoment])) #"odontophoreangle"
    ocenter = float(position[discretemoment]) #"x"
    
    topcontactpointy = float(ytop[discretemoment])
    bottomcontactpointy = float(ybottom[discretemoment])
    # bottomofodontophorex = (ocenter - (oradius*np.cos(np.deg2rad(float(angle[discretemoment])))))
    # bottomofodontophorey = (0 - (oradius*np.sin(np.deg2rad(float(angle[discretemoment])))))
    #sideofodontophorex = (ocenter - (ominoraxisradius*np.cos(np.deg2rad(float(angle[discretemoment])))))
    #sideofodontophorey = 0
    
    i1i3radius = float(i1i3radiusarray[discretemoment])
    i2length = float(i2lengtharray[discretemoment])
    topangle = float(topanglearray[discretemoment]) #Tphi
    bottomangle = float(bottomanglearray[discretemoment]) #Bphi
    xtopval = -(np.cos(np.deg2rad(topangle))*i2length)/2
    xbottomval = -(np.cos(np.deg2rad(bottomangle))*i2length)/2
    
    furthestbackxpoint = float(furthestbackxpointarray[discretemoment])
    furthestbackypoint = float(furthestbackypointarray[discretemoment])
    i1i3contacttopy = float(i1i3contacttopyarray[discretemoment])
    i1i3contacttopx = -0.00125
    i1i3contactbottomy = float(i1i3contactbottomyarray[discretemoment])
    i1i3contactbottomx = -.00125

    #x1val = float(x1array[discretemoment])

    ocontacttopx = float(ocontacttopxarray[discretemoment])
    ocontacttopy = float(ocontacttopyarray[discretemoment])
    ocontactbottomx = float(ocontactbottomxarray[discretemoment])
    ocontactbottomy = float(ocontactbottomyarray[discretemoment])

    #bigxval = float(bigxvalarray[discretemoment])
    
    #Create Shapes
    #odontophore
    odontophore = Ellipse((ocenter, 0), ominoraxisdiameter, odiameter, oangle)
    odontophore.set_alpha(1)
    odontophore.set_edgecolor('black')
    odontophore.set_facecolor('none')
    odontophore.set_lw(1)
    
    #i1i3
    i1i3top = Circle((0,topcontactpointy),0.00125)
    i1i3top.set_edgecolor('black')
    i1i3top.set_facecolor('none')
    i1i3bottom = Circle((0,bottomcontactpointy),0.00125)
    i1i3bottom.set_edgecolor('black')
    i1i3bottom.set_facecolor('none')
    
    #i2
    #Top line
    x1,y1 = np.array([[i1i3contacttopx, ocontacttopx],[i1i3contacttopy, ocontacttopy]])
    #Bottome line
    x2, y2 = np.array([[i1i3contactbottomx, ocontactbottomx],[i1i3contactbottomy, ocontactbottomy]])
    i2top = mlines.Line2D(x1,y1)
    i2bottom= mlines.Line2D(x2,y2)
    
    #Add each shape to the figure
    figure.add_artist(odontophore)
    figure.add_artist(i1i3top)
    figure.add_artist(i1i3bottom)
    figure.add_artist(i2top)
    figure.add_artist(i2bottom)
    plt.xlim(-.02,.02)
    plt.ylim(-.01,.01)
    return figure

#Important function to erase shapes from the figure, otherwise each frame would overlap over eachother
def resetVar(a):
    a.remove()
    a = plt.subplot(111, aspect='equal')
    return a


######################################################
##                 Create Video                     ##
######################################################


#set which movie writer class is being used
if animation.FFMpegWriter.isAvailable():
    movieWriterClass = animation.FFMpegWriter
#elif animation.AVConvWvriter.isAvailable():
    #movieWriterClass = animation.AVConvWriter
else:
    sys.stderr.write('video converter missing')
    exit()

#For each discretemoment/time step, save the figure as a frame in the .mp4 file
moviewriter = movieWriterClass(fps = fps, bitrate = br)
with moviewriter.saving(plt.gcf(), moviefile, dpi):
    i = 1 #first frame, must be 1<i<850
    frames = 850 #must be frames<=850
    while( i < frames):
        sys.stdout.flush()
        progress = str(round((i*100)/frames))
        sys.stdout.write(progress + '% complete \r' )
        createShapes(i, a)
        moviewriter.grab_frame()
        a = resetVar(a)
        i = i+1
moviewriter.finish()

