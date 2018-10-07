import os
import csv
import matplotlib
matplotlib.use('Agg') # force matplotlib to not use XWindows backend
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import animation
#import plotly.plotly as py
#import plotly.tools as tls

#The following code reads the SlugOutput2.csv and saves the time in an array list "time"
time = []
i = 0
with open("SlugOutput2.csv", "r") as file:
     row = csv.reader(file)
     while i < 850:
        for column in row:
            time.insert(i,column[0])
            i+=1
time[1:] = list(map(float, time[1:]))

#The following line prints the time array to confirm that it contains the correct values
#print(time)

#The following code reads the SlugOutput2.csv and saves the position in an array list "position"
position = []
i = 0
with open("SlugOutput2.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            position.insert(i,column[1])
            i+=1
position[1:] = list(map(float, position[1:]))

#The following line prints the time array to confirm that it contains the correct values
#print(position)

#The following code reads the SlugOutput2.csv and saves the radius in an array list "radius"
radius = []
i = 0
with open("SlugOutput2.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            radius.insert(i,column[2])
            i+=1
radius[1:] = list(map(float, radius[1:]))

#The following line prints the time array to confirm that it contains the correct values
#print(radius)

#The following code reads the SlugOutput2.csv and saves the angle in an array list "angle"
angle = []
i = 0
with open("SlugOutput2.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            angle.insert(i,column[3])
            i+=1
angle[1:] = list(map(float, angle[1:]))

#The following line prints the time array to confirm that it contains the correct values
#print(angle)

#The following code reads the SlugOutput2.csv and saves the seaweedforce in an array list "seaweedforce"
seaweedforce = []
i = 0
with open("SlugOutput2.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            seaweedforce.insert(i,column[14])
            i+=1
seaweedforce[1:] = list(map(float, seaweedforce[1:]))

#The following line prints the time array to confirm that it contains the correct values
#print(seaweedforce)

#The following code reads the SlugOutput2.csv and saves the Odontophore input in an array list "freqN3"
freqN3 = []
i = 0
with open("SlugOutput2.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            freqN3.insert(i,column[8])
            i+=1
freqN3[1:] = list(map(float, freqN3[1:]))

#The following line prints the time array to confirm that it contains the correct values
#print(freqN3)

#The following code reads the SlugOutput2.csv and saves the Hinge input in an array list "freqHinge"
freqHinge = []
i = 0
with open("SlugOutput2.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            freqHinge.insert(i,column[9])
            i+=1
freqHinge[1:] = list(map(float, freqHinge[1:]))

#The following line prints the time array to confirm that it contains the correct values
#print(freqHinge)

#The following code reads the SlugOutput2.csv and saves the I2 input in an array list "freqI2"
freqI2 = []
i = 0
with open("SlugOutput2.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            freqI2.insert(i,column[6])
            i+=1
freqI2[1:] = list(map(float, freqI2[1:]))

#The following line prints the time array to confirm that it contains the correct values
#print(freqI2)

#The following code reads the SlugOutput2.csv and saves the I1/I3 input in an array list "freqI1I3"
freqI1I3 = []
i = 0
with open("SlugOutput2.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            freqI1I3.insert(i,column[7])
            i+=1
freqI1I3[1:] = list(map(float, freqI1I3[1:]))

#The following line prints the time array to confirm that it contains the correct values
#print(freqI1I3)

#The following code plots the entire figure
figure = plt.figure(figsize = (18,18))
timeticks = np.arange(0,9,1)
#The figure is a Grid with 5 rows, 2 columns. The input graphs are in the left column, output in the right
gs = gridspec.GridSpec(5, 2)

#The following code subplots the position vs time graph
positiongraph = figure.add_subplot(gs[0,1])
positionticks = np.arange(-.006,.004, .001)
positiongraph.plot(time[1:850],position[1:850])
positiongraph.set_title('Odontophore Position')
positiongraph.set_xlabel(time[0] + ' (seconds)')
positiongraph.set_ylabel(position[0] + ' (Meters)')
positiongraph.axis([0, 9, -0.006, 0.004])
positiongraph.set_xticks(timeticks)
positiongraph.set_yticks(positionticks)
positiongraph.grid(True)

#The following code subplots the Odontophore major axis vs time graph
radiusgraph = figure.add_subplot(gs[1,1])
radiusticks = np.arange(0,.009, .001)
radiusgraph.plot(time[1:850],radius[1:850])
radiusgraph.set_title('Odontophore Major Axis')
radiusgraph.set_xlabel(time[0] + ' (seconds)')
radiusgraph.set_ylabel(radius[0] + ' (Meters)')
radiusgraph.axis([0, 9, 0, 0.009])
radiusgraph.set_xticks(timeticks)
radiusgraph.set_yticks(radiusticks)
radiusgraph.grid(True)

#The following code subplots the Angle vs time graph
anglegraph = figure.add_subplot(gs[2,1])
angleticks = np.arange(-10,100, 10)
anglegraph.plot(time[1:850],angle[1:850])
anglegraph.set_title('Odontophore Angle')
anglegraph.set_xlabel(time[0] + ' (seconds)')
anglegraph.set_ylabel(angle[0] + ' (Degrees)')
anglegraph.axis([0, 9, -10, 100])
anglegraph.set_xticks(timeticks)
anglegraph.set_yticks(angleticks)
anglegraph.grid(True)

#The following code subplots the Seaweed Force vs time graph
forcegraph = figure.add_subplot(gs[0,0])
forceticks = np.arange(-20,120, 20)
forcegraph.plot(time[1:850],seaweedforce[1:850])
forcegraph.set_title('Force From Seaweed')
forcegraph.set_xlabel(time[0] + ' (seconds)')
forcegraph.set_ylabel(seaweedforce[0] + ' (mN)')
forcegraph.axis([0, 9, -20, 120])
forcegraph.set_xticks(timeticks)
forcegraph.set_yticks(forceticks)
forcegraph.grid(True)

#The following code subplots the Odontophore input vs time graph
n3graph = figure.add_subplot(gs[1,0])
n3ticks = np.arange(-5,35, 5)
n3graph.plot(time[1:850],freqN3[1:850])
n3graph.set_title('Odontophore Input')
n3graph.set_xlabel(time[0] + ' (seconds)')
n3graph.set_ylabel(freqN3[0] + ' (Hz)')
n3graph.axis([0, 9, -5, 35])
n3graph.set_xticks(timeticks)
n3graph.set_yticks(n3ticks)
n3graph.grid(True)

#The following code subplots the Hinge input vs time graph
hingegraph = figure.add_subplot(gs[2,0])
hingeticks = np.arange(-5,25, 5)
hingegraph.plot(time[1:850],freqHinge[1:850])
hingegraph.set_title('Hinge Input')
hingegraph.set_xlabel(time[0] + ' (seconds)')
hingegraph.set_ylabel(freqHinge[0] + ' (Hz)')
hingegraph.axis([0, 9, -5, 25])
hingegraph.set_xticks(timeticks)
hingegraph.set_yticks(hingeticks)
hingegraph.grid(True)

#The following code subplots the I2 input vs time graph
i2graph = figure.add_subplot(gs[3,0])
i2ticks = np.arange(-5,25,5)
i2graph.plot(time[1:850],freqI2[1:850])
i2graph.set_title('I2 Input')
i2graph.set_xlabel(time[0] + ' (seconds)')
i2graph.set_ylabel(freqI2[0] + ' (Hz)')
i2graph.axis([0, 9, -5, 25])
i2graph.set_xticks(timeticks)
i2graph.set_yticks(i2ticks)
i2graph.grid(True)

#The following code subplots the I1/I3 vs time graph
i1i3graph = figure.add_subplot(gs[4,0])
i1i3ticks = np.arange(-5,25,5)
i1i3graph.plot(time[1:850],freqI1I3[1:850])
i1i3graph.set_title('I1/I3 Input')
i1i3graph.set_xlabel(time[0] + ' (seconds)')
i1i3graph.set_ylabel(freqI1I3[0] + ' (Hz)')
i1i3graph.axis([0, 9, -5, 25])
i1i3graph.set_xticks(timeticks)
i1i3graph.set_yticks(i1i3ticks)
i1i3graph.grid(True)

#Show the entire figure
#plt.show()

#Save Figure to PDF
figure.tight_layout()
figure.savefig("plot.pdf")



#Design the 2D animation of the Buccal Mass

#first an example


# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))
line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = np.linspace(0, 2, 1000)
    y = np.sin(2 * np.pi * (x - 0.01 * i))
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=200, interval=20, blit=True)

#plt.show()
