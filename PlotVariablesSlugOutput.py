import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg') # force matplotlib to not use XWindows backend
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# read in the data
with open("SlugOutput2.csv","r") as file:
    data = pd.read_csv(file)

#The following code plots the entire figure
figure = plt.figure(figsize = (18,18))
timeticks = np.arange(0,9,1)

#The figure is a Grid with 5 rows, 2 columns. The input graphs are in the left column, output in the right
gs = gridspec.GridSpec(5, 2)

#The following code subplots the position vs time graph
positiongraph = figure.add_subplot(gs[0,1])
positionticks = np.arange(-.006,.004, .001)
positiongraph.plot(data['time'], data['position'])
positiongraph.set_title('Odontophore Position')
positiongraph.set_xlabel('Time (seconds)')
positiongraph.set_ylabel('Position (meters)')
positiongraph.axis([0, 9, -0.006, 0.004])
positiongraph.set_xticks(timeticks)
positiongraph.set_yticks(positionticks)
positiongraph.grid(True)

#The following code subplots the Odontophore major axis vs time graph
radiusgraph = figure.add_subplot(gs[1,1])
radiusticks = np.arange(0,.009, .001)
radiusgraph.plot(data['time'], data['radius'])
radiusgraph.set_title('Odontophore Major Axis')
radiusgraph.set_xlabel('Time (seconds)')
radiusgraph.set_ylabel('Radius (meters)')
radiusgraph.axis([0, 9, 0, 0.009])
radiusgraph.set_xticks(timeticks)
radiusgraph.set_yticks(radiusticks)
radiusgraph.grid(True)

#The following code subplots the Angle vs time graph
anglegraph = figure.add_subplot(gs[2,1])
angleticks = np.arange(-10,100, 10)
anglegraph.plot(data['time'], data['angle'])
anglegraph.set_title('Odontophore Angle')
anglegraph.set_xlabel('Time (seconds)')
anglegraph.set_ylabel('Angle (degrees)')
anglegraph.axis([0, 9, -10, 100])
anglegraph.set_xticks(timeticks)
anglegraph.set_yticks(angleticks)
anglegraph.grid(True)

#The following code subplots the Seaweed Force vs time graph
forcegraph = figure.add_subplot(gs[0,0])
forceticks = np.arange(-20,120, 20)
forcegraph.plot(data['time'], data['seaweedforce'])
forcegraph.set_title('Force From Seaweed')
forcegraph.set_xlabel('Time (seconds)')
forcegraph.set_ylabel('Seaweed Force (mN)')
forcegraph.axis([0, 9, -20, 120])
forcegraph.set_xticks(timeticks)
forcegraph.set_yticks(forceticks)
forcegraph.grid(True)

#The following code subplots the Odontophore input vs time graph
n3graph = figure.add_subplot(gs[1,0])
n3ticks = np.arange(-5,35, 5)
n3graph.plot(data['time'], data['freqN3'])
n3graph.set_title('Odontophore Input')
n3graph.set_xlabel('Time (seconds)')
n3graph.set_ylabel('Odontophore Frequency (Hz)')
n3graph.axis([0, 9, -5, 35])
n3graph.set_xticks(timeticks)
n3graph.set_yticks(n3ticks)
n3graph.grid(True)

#The following code subplots the Hinge input vs time graph
hingegraph = figure.add_subplot(gs[2,0])
hingeticks = np.arange(-5,25, 5)
hingegraph.plot(data['time'], data['freqHinge'])
hingegraph.set_title('Hinge Input')
hingegraph.set_xlabel('Time (seconds)')
hingegraph.set_ylabel('Hinge Frequency (Hz)')
hingegraph.axis([0, 9, -5, 25])
hingegraph.set_xticks(timeticks)
hingegraph.set_yticks(hingeticks)
hingegraph.grid(True)

#The following code subplots the I2 input vs time graph
i2graph = figure.add_subplot(gs[3,0])
i2ticks = np.arange(-5,25,5)
i2graph.plot(data['time'], data['freqI2'])
i2graph.set_title('I2 Input')
i2graph.set_xlabel('Time (seconds)')
i2graph.set_ylabel('I2 Frequency (Hz)')
i2graph.axis([0, 9, -5, 25])
i2graph.set_xticks(timeticks)
i2graph.set_yticks(i2ticks)
i2graph.grid(True)

#The following code subplots the I1/I3 vs time graph
i1i3graph = figure.add_subplot(gs[4,0])
i1i3ticks = np.arange(-5,25,5)
i1i3graph.plot(data['time'], data['freqI1I3'])
i1i3graph.set_title('I1/I3 Input')
i1i3graph.set_xlabel('Time (seconds)')
i1i3graph.set_ylabel('I1/I3 Frequency (Hz)')
i1i3graph.axis([0, 9, -5, 25])
i1i3graph.set_xticks(timeticks)
i1i3graph.set_yticks(i1i3ticks)
i1i3graph.grid(True)

#Show the entire figure
#plt.show()

#Save Figure to PDF
figure.tight_layout()
figure.savefig("plot.pdf")
