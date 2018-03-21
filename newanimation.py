import os
import sys
import csv
import matplotlib
#matplotlib.use('Agg') # force matplotlib to not use XWindows backend
from math import sqrt
import matplotlib.gridspec as gridspec
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
height = 800
#frames per second
fps = 25
#bitrate to reduce compression
br = 2000
#moviesize
dpi = 180


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

#The following code reads the animationinfo-csvTranslated.csv and saves the time in an array list "time"
time = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            time.insert(i,column[0])
            i+=1

#The following line prints the time array to confirm that it contains the correct values
#print(time)

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
            bigxvalarray.insert(i,column[20])
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
            x1array.insert(i,column[21])
            i+=1

#The following line prints the x1array to confirm that it contains the correct values
#print(x1array)

#The following code reads the animationinfo-csvTranslated.csv and saves the freqi2 in an array list "freqi2array"
freqi2array = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            freqi2array.insert(i,column[22])
            i+=1

#The following line prints the freqi2array to confirm that it contains the correct values
#print(freqi2array)

#The following code reads the animationinfo-csvTranslated.csv and saves the freqi1i3 in an array list "freqi1i3array"
freqi1i3array = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            freqi1i3array.insert(i,column[23])
            i+=1

#The following line prints the freqi1i3array to confirm that it contains the correct values
#print(freqi1i3array)

#The following code reads the animationinfo-csvTranslated.csv and saves the freqN3 in an array list "freqN3array"
freqN3array = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            freqN3array.insert(i,column[24])
            i+=1

#The following line prints the freqN3array to confirm that it contains the correct values
#print(freqN3array)

#The following code reads the animationinfo-csvTranslated.csv and saves the freqHinge in an array list "freqHingearray"
freqHingearray = []
i = 0
with open("animationinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            freqHingearray.insert(i,column[25])
            i+=1

#The following line prints the freqHingearray to confirm that it contains the correct values
#print(freqHingearray)

######################################################
##                Plot 2D Simulation                ##
######################################################


figure = plt.figure()
gs = gridspec.GridSpec(18, 9)

# graphs for frequencies, a is animation figure
i2graph = figure.add_subplot(gs[11:12,1:8])
n3graph = figure.add_subplot(gs[13:14,1:8])
i1i3graph = figure.add_subplot(gs[15:16,1:8])
hingegraph = figure.add_subplot(gs[17:18,1:8])
a = figure.add_subplot(gs[2:9,0:9], aspect='equal')

# shapes
odontophore = Ellipse((0, 0), 0, 0, 0)
i1i3top = Circle((0,0),0.00125)
i1i3bottom = Circle((0,0),0.00125)
i2top = mlines.Line2D([0,1],[0,1])
i2bottom= mlines.Line2D([0,1],[0,1])
    
a.add_artist(odontophore)
a.add_artist(i1i3top)
a.add_artist(i1i3bottom)
a.add_artist(i2top)
a.add_artist(i2bottom)
plt.ylabel('Y- Plane (m)', fontsize = 4)
plt.xlabel('X- Plane (m)', fontsize = 4)
plt.xlim(-.02,.02)
plt.ylim(-.01,.01)
plt.title('Animation', fontsize = 8)
a.tick_params(labelsize=3)
a.text(-.0240,.0175, 'Channels Activated:', color='Black', fontsize=4)
a.text(-.0240,.015, 'Level of Activation (Hz):', color='Black', fontsize=4)

#a text
toptextI2 = a.text(.0140,.0175, '- I 2 -', color='black', fontsize=4)
toptextN3 = a.text(-.008,.0175, '- N 3 -', color='black', fontsize=4)
toptexti1i3 = a.text(-.0020,.0175, '- I 1 I 3 -', color='black', fontsize=4)
toptexthinge = a.text(.0045,.0175, '- H I N G E -', color='black', fontsize=4)
bottomtexti2 = a.text(.0140,.015, '%2.1f' % float(freqi2array[1]), color='black', fontsize=4)
bottomtextN3 = a.text(-.008,.015, '%2.1f' % float(freqN3array[1]), color='black', fontsize=4)
bottomtexti1i3 = a.text(-.0020,.015, '%2.1f' % float(freqi1i3array[1]), color='black', fontsize=4)
bottomtexthinge = a.text(.0045,.015, '%2.1f' % float(freqHingearray[1]), color='black', fontsize=4)
timer = a.text(0, .02, 'Time: 0 seconds', color='Black', fontsize=6, horizontalalignment = 'center')

timeticks = np.arange(0,9,1)
#The following code subplots the I2 input vs time graph
i2ticks = np.arange(-5,35,5)
i2graph.plot(time[1:850],freqi2array[1:850])
#i2graph.set_title('I2 Input', fontsize = 4)
#i2graph.set_xlabel(time[0] + ' (seconds)', fontsize = 2)
i2graph.set_ylabel(freqi2array[0] + ' (Hz)', fontsize = 2)
i2graph.axis([0, 9, -5, 35])
i2graph.set_xticks(timeticks)
i2graph.set_yticks(i2ticks)
i2graph.tick_params(labelsize=3)
i2graph.set_axis_off()
i2graph.text(9.0140,10, '- I 2 -', color='black', fontsize=4)
i2graph.text(-0.5140,30, '30', color='black', fontsize=4)
i2graph.text(-0.5140,00, '0', color='black', fontsize=4)
#i2graph.grid(True)
    
#The following code subplots the Odontophore input vs time graph
n3ticks = np.arange(-5,40, 5)
n3graph.plot(time[1:850],freqN3array[1:850])
#n3graph.set_title('Odontophore Input', fontsize = 4)
#n3graph.set_xlabel(time[0] + ' (seconds)', fontsize = 2)
n3graph.set_ylabel(freqN3array[0] + ' (Hz)', fontsize = 2)
n3graph.axis([0, 9, -5, 40], fontsize = 2)
n3graph.set_xticks(timeticks)
n3graph.set_yticks(n3ticks)
n3graph.tick_params(labelsize=3)
n3graph.set_axis_off()
n3graph.text(9.0140,10, '- N 3 -', color='black', fontsize=4)
n3graph.text(-0.5140,35, '35', color='black', fontsize=4)
n3graph.text(-0.5140,00, '0', color='black', fontsize=4)
#n3graph.grid(True)
#The following code subplots the I1/I3 vs time graph
i1i3ticks = np.arange(-5,75,5)
i1i3graph.plot(time[1:850],freqi1i3array[1:850])
#i1i3graph.set_title('I1/I3 Input', fontsize = 4)
i1i3graph.set_xlabel(time[0] + ' (seconds)', fontsize = 2)
i1i3graph.set_ylabel(freqi1i3array[0] + ' (Hz)', fontsize = 2)
i1i3graph.axis([0, 9, -5, 75])
i1i3graph.set_xticks(timeticks)
i1i3graph.set_yticks(i1i3ticks)
i1i3graph.tick_params(labelsize=3)
i1i3graph.set_axis_off()
i1i3graph.text(9.0140,10, '- I 1 I 3 -', color='black', fontsize=4)
i1i3graph.text(-0.5140,65, '65', color='black', fontsize=4)
i1i3graph.text(-0.5140,00, '0', color='black', fontsize=4)
#i1i3graph.grid(True)

#The following code subplots the Hinge input vs time graph
hingeticks = np.arange(-5,15, 5)
hingegraph.plot(time[1:850],freqHingearray[1:850])
#hingegraph.set_title('Hinge Input', fontsize = 4)
#hingegraph.set_xlabel(time[0] + ' (seconds)', fontsize = 2)
hingegraph.set_ylabel(freqHingearray[0] + ' (Hz)', fontsize = 2)
hingegraph.axis([0, 9, -5, 15], fontsize = 2)
hingegraph.set_xticks(timeticks)
hingegraph.set_yticks(hingeticks)
hingegraph.tick_params(labelsize=3)
hingegraph.set_axis_off()
hingegraph.text(9.0140,10, '- HINGE -', color='black', fontsize=4)
hingegraph.text(-0.5140,10, '10', color='black', fontsize=4)
hingegraph.text(-0.5140,00, '0', color='black', fontsize=4)
#hingegraph.grid(True)

vertlinei2 = i2graph.vlines([float(time[1])],[0],80, colors = 'r')
vertlinen3 = n3graph.vlines([float(time[1])],[0],80, colors = 'r')
vertlinei1i3 = i1i3graph.vlines([float(time[1])],[0],80, colors = 'r')
vertlinehinge = hingegraph.vlines([float(time[1])],[0],80, colors = 'r')

def makevertlinei2(i2graph, discretemoment):
    return i2graph.vlines([float(time[discretemoment])],[0],80, colors = 'r')
    
def makevertlinen3(n3graph, discretemoment):
    return n3graph.vlines([float(time[discretemoment])],[0],80, colors = 'r')
    
def makevertlinei1i3(i1i3graph, discretemoment):
    return i1i3graph.vlines([float(time[discretemoment])],[0],80, colors = 'r')
    
def makevertlinehinge(hingegraph, discretemoment):
    return hingegraph.vlines([float(time[discretemoment])],[0],80, colors = 'r')

def createOdontophore(discretemoment, odontophore):
    #Define Variables by acsessing array values
    oradius = float(radius[discretemoment]) #"a" in the kinetic model
    odiameter = 2*oradius #2*a
    ominoraxisradius = (5.0 * sqrt(5.0) / sqrt(oradius * 1000.0)) / 1000.0 #"b" in the kinetic model
    ominoraxisdiameter = 2*ominoraxisradius #2*b
    oangle = -(90 - float(angle[discretemoment])) #"odontophoreangle"
    ocenter = float(position[discretemoment]) #"x"
    #odontophore
    odontophore = Ellipse((ocenter, 0), ominoraxisdiameter, odiameter, oangle)
    odontophore.set_alpha(1)
    odontophore.set_facecolor('none')
    odontophore.set_linewidth(.5)
    if float(freqN3array[discretemoment]) == 0:
        odontophore.set_edgecolor('black')
    elif float(freqN3array[discretemoment]) > 0:
        odontophore.set_edgecolor('black')
    else:
        odontophore.set_edgecolor('red')
    return odontophore

def createi1i3top(discretemoment, i1i3top):
    topcontactpointy = float(ytop[discretemoment])
    i1i3top = Circle((0,topcontactpointy),0.00125)
    i1i3top.set_facecolor('none')
    i1i3top.set_linewidth(.5)
    if float(freqi1i3array[discretemoment]) == 0:
        i1i3top.set_edgecolor('black')
    elif float(freqi1i3array[discretemoment]) > 0:
        i1i3top.set_edgecolor('black')
    else:
        i1i3top.set_edgecolor('red')
    return i1i3top

def createi1i3bottom(discretemoment, i1i3bottom):
    bottomcontactpointy = float(ybottom[discretemoment])
    i1i3bottom = Circle((0,bottomcontactpointy),0.00125)
    i1i3bottom.set_facecolor('none')
    i1i3bottom.set_linewidth(.5)
    if float(freqi1i3array[discretemoment]) == 0:
        i1i3bottom.set_edgecolor('black')
    elif float(freqi1i3array[discretemoment]) > 0:
        i1i3bottom.set_edgecolor('black')
    else:
        i1i3bottom.set_edgecolor('red')
    return i1i3bottom

def createi2top(discretemoment, i2top):
    #Top line
    
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

    x1,y1 = np.array([[i1i3contacttopx, ocontacttopx],[i1i3contacttopy, ocontacttopy]])
    i2top = mlines.Line2D(x1,y1)
    i2top.set_linewidth(.5)
    if float(freqi2array[discretemoment]) == 0:
        i2top.set_color('black')
    elif float(freqi2array[discretemoment]) > 0:
        i2top.set_color('black')
    else:
        i2top.set_color('red')
    return i2top

def createi2bottom(discretemoment, i2bottom):
    #Bottom line
    
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

    x2, y2 = np.array([[i1i3contactbottomx, ocontactbottomx],[i1i3contactbottomy, ocontactbottomy]])
    i2bottom= mlines.Line2D(x2,y2)
    i2bottom.set_linewidth(.5)
    if float(freqi2array[discretemoment]) == 0:
        i2bottom.set_color('black')
    elif float(freqi2array[discretemoment]) > 0:
        i2bottom.set_color('black')
    else:
        i2bottom.set_color('red')
    return i2bottom

def maketoptextI2(figure, discretemoment, toptexti2):
    if float(freqi2array[discretemoment]) == 0:
        toptexti2 = figure.text(.0140,.0175, '- I 2 -', color='black', fontsize=4)
    elif float(freqi2array[discretemoment]) > 0:
        toptexti2 = figure.text(.0140,.0175, '- I 2 -', color='green', fontsize=4)
    else:
        toptexti2 = figure.text(.0140,.0175, '%f' % float(freqi2array[discretemoment]), color='red', fontsize=4)
    return toptexti2

def maketoptextN3(figure, discretemoment, toptextN3):
    if float(freqN3array[discretemoment]) == 0:
        toptextN3 = figure.text(-.008,.0175, '- N 3 -', color='black', fontsize=4)
    elif float(freqN3array[discretemoment]) > 0:
        toptextN3 = figure.text(-.008,.0175, '- N 3 -', color='green', fontsize=4)
    else:
        toptextN3 = figure.text(-.008,.0175, '%f' % float(freqN3array[discretemoment]), color='red', fontsize=4)
    return toptextN3

def maketoptexti1i3(figure, discretemoment, toptexti1i3):
    if float(freqi1i3array[discretemoment]) == 0:
        toptexti1i3 = figure.text(-.0020,.0175, '- I 1 I 3 -', color='black', fontsize=4)
    elif float(freqi1i3array[discretemoment]) > 0:
        toptexti1i3 = figure.text(-.0020,.0175, '- I 1 I 3 -', color='green', fontsize=4)
    else:
        toptexti1i3 = figure.text(-.0020,.0175, '%f' % float(freqi1i3array[discretemoment]), color='red', fontsize=4)
    return toptexti1i3

def maketoptexthinge(figure, discretemoment, toptexthinge):
    if float(freqHingearray[discretemoment]) == 0:
        toptexthinge = figure.text(.0045,.0175, '- H I N G E -', color='black', fontsize=4)
    elif float(freqHingearray[discretemoment]) > 0:
        toptexthinge = figure.text(.0045,.0175, '- H I N G E -', color='green', fontsize=4)
    else:
        toptexthinge = figure.text(.0025,.0175, '%f' % float(freqHingearray[discretemoment]), color='red', fontsize=4)
    return toptexthinge

def makebottomtexti2(figure, discretemoment, bottomtexti2):
    if float(freqi2array[discretemoment]) == 0:
        bottomtexti2 = figure.text(.0140,.015, '%2.0f' % float(freqi2array[discretemoment]), color='black', fontsize=4)
    elif float(freqi2array[discretemoment]) > 0:
        bottomtexti2 = figure.text(.0140,.015, '%2.0f' % float(freqi2array[discretemoment]), color='green', fontsize=4)
    else:
        bottomtexti2 = figure.text(.0140,.015, '%f' % float(freqi2array[discretemoment]), color='red', fontsize=4)
    return bottomtexti2

def makebottomtextN3(figure, discretemoment, bottomtextN3):
    if float(freqN3array[discretemoment]) == 0:
        bottomtextN3 = figure.text(-.008,.015, '%2.0f' % float(freqN3array[discretemoment]), color='black', fontsize=4)
    elif float(freqN3array[discretemoment]) > 0:
        bottomtextN3 = figure.text(-.008,.015, '%2.0f' % float(freqN3array[discretemoment]), color='green', fontsize=4)
    else:
        bottomtextN3 = figure.text(-.008,.015, '%f' % float(freqN3array[discretemoment]), color='red', fontsize=4)
    return bottomtextN3

def makebottomtexti1i3(figure, discretemoment, bottomtexti1i3):
    if float(freqi1i3array[discretemoment]) == 0:
        bottomtexti1i3 = figure.text(-.0020,.015, '%2.0f' % float(freqi1i3array[discretemoment]), color='black', fontsize=4)
    elif float(freqi1i3array[discretemoment]) > 0:
        bottomtexti1i3 = figure.text(-.0020,.015, '%2.0f' % float(freqi1i3array[discretemoment]), color='green', fontsize=4)
    else:
        bottomtexti1i3 = figure.text(-.0020,.015, '%f' % float(freqi1i3array[discretemoment]), color='red', fontsize=4)
    return bottomtexti1i3

def makebottomtexthinge(figure, discretemoment, bottomtexthinge):
    if float(freqHingearray[discretemoment]) == 0:
        bottomtexthinge = figure.text(.0045,.015, '%2.0f' % float(freqHingearray[discretemoment]), color='black', fontsize=4)
    elif float(freqHingearray[discretemoment]) > 0:
        bottomtexthinge = figure.text(.0045,.015, '%2.0f' % float(freqHingearray[discretemoment]), color='green', fontsize=4)
    else:
        bottomtexthinge = figure.text(.0045,.015, '%f' % float(freqHingearray[discretemoment]), color='red', fontsize=4)
    return bottomtexthinge

def maketimer(figure, discretemoment, timer):
    timer = figure.text(0, .02, 'Time: %3.2f seconds' % float(discretemoment/100.00), color='Black', fontsize=6, horizontalalignment = 'center')
    return timer

def resetn3(n3graph):
    n3graph.remove()
    n3graph = figure.add_subplot(gs[13:14,1:8])
    #plt.tight_layout()
    return n3graph

def resethinge(hingegraph):
    hingegraph.remove()
    hingegraph = figure.add_subplot(gs[17:18,1:8])
    #plt.tight_layout()
    return hingegraph

def reseti2(i2graph):
    i2graph.remove()
    i2graph = figure.add_subplot(gs[11:12,1:8])
    #plt.tight_layout()
    return i2graph

def reseti1i3(i1i3graph):
    i1i3graph.remove()
    i1i3graph = figure.add_subplot(gs[15:16,1:8])
    #plt.tight_layout()
    return i1i3graph

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
        #makefigures(i2graph, n3graph, i1i3graph, hingegraph, i)
        odontophore.remove()
        i1i3top.remove()
        i1i3bottom.remove()
        i2top.remove()
        i2bottom.remove()
        vertlinei2.remove()
        vertlinen3.remove()
        vertlinei1i3.remove()
        vertlinehinge.remove()
        toptextI2.remove()
        toptextN3.remove()
        toptexti1i3.remove()
        toptexthinge.remove()
        bottomtexti2.remove()
        bottomtextN3.remove()
        bottomtexti1i3.remove()
        bottomtexthinge.remove()
        timer.remove()
        toptextI2 = maketoptextI2(a, i, toptextI2)
        toptextN3 = maketoptextN3(a, i, toptextN3)
        toptexti1i3 = maketoptexti1i3(a, i, toptexti1i3)
        toptexthinge = maketoptexthinge(a, i, toptexthinge)
        bottomtexti2 = makebottomtexti2(a, i, bottomtexti2)
        bottomtextN3 = makebottomtextN3(a, i, bottomtextN3)
        bottomtexti1i3 = makebottomtexti1i3(a, i, bottomtexti1i3)
        bottomtexthinge = makebottomtexthinge(a, i, bottomtexthinge)
        timer = maketimer(a, i, timer)
        vertlinei2 = makevertlinei2(i2graph, i)
        vertlinen3 = makevertlinen3(n3graph, i)
        vertlinei1i3 = makevertlinei1i3(i1i3graph, i)
        vertlinehinge = makevertlinehinge(hingegraph, i)
        odontophore = createOdontophore(i, odontophore)
        i1i3top = createi1i3top(i, i1i3top)
        i1i3bottom = createi1i3bottom(i, i1i3bottom)
        i2top = createi2top(i, i2top)
        i2bottom = createi2bottom(i, i2bottom)
        a.add_artist(odontophore)
        a.add_artist(i1i3top)
        a.add_artist(i1i3bottom)
        a.add_artist(i2top)
        a.add_artist(i2bottom)
        moviewriter.grab_frame()
        i = i+4
moviewriter.finish()
