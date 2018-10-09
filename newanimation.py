import sys
from math import sqrt
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg') # force matplotlib to not use XWindows backend
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
from matplotlib.patches import Ellipse
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
fps = 20
#bitrate to reduce compression
br = 2000
#moviesize
dpi = 180


######################################################
##                    Read File                     ##
######################################################

# read in the data
with open("animationinfo.csv","r") as file:
    data = pd.read_csv(file)

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
toptextN3 = a.text(-.008,.0175, '- I 4 -', color='black', fontsize=4)
toptexti1i3 = a.text(-.0020,.0175, '- I 1 I 3 -', color='black', fontsize=4)
toptexthinge = a.text(.0045,.0175, '- H I N G E -', color='black', fontsize=4)
bottomtexti2 = a.text(.0140,.015, '%2.1f' % data['freqI2'][0], color='black', fontsize=4)
bottomtextN3 = a.text(-.008,.015, '%2.1f' % data['freqN3'][0], color='black', fontsize=4)
bottomtexti1i3 = a.text(-.0020,.015, '%2.1f' % data['freqI1I3'][0], color='black', fontsize=4)
bottomtexthinge = a.text(.0045,.015, '%2.1f' % data['freqHinge'][0], color='black', fontsize=4)
timer = a.text(0, .02, 'Time: 0 seconds', color='Black', fontsize=6, horizontalalignment = 'center')

timeticks = np.arange(0,9,1)
#The following code subplots the I2 input vs time graph
i2ticks = np.arange(-5,35,5)
i2graph.plot(data['time'], data['freqI2'])
#i2graph.set_title('I2 Input', fontsize = 4)
#i2graph.set_xlabel('Time (seconds)', fontsize = 2)
i2graph.set_ylabel('I2 Frequency (Hz)', fontsize = 2)
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
n3graph.plot(data['time'], data['freqN3'])
#n3graph.set_title('Odontophore Input', fontsize = 4)
#n3graph.set_xlabel('Time (seconds)', fontsize = 2)
n3graph.set_ylabel('Odontophore Frequency (Hz)', fontsize = 2)
n3graph.axis([0, 9, -5, 40], fontsize = 2)
n3graph.set_xticks(timeticks)
n3graph.set_yticks(n3ticks)
n3graph.tick_params(labelsize=3)
n3graph.set_axis_off()
n3graph.text(9.0140,10, '- I 4 -', color='black', fontsize=4)
n3graph.text(-0.5140,35, '35', color='black', fontsize=4)
n3graph.text(-0.5140,00, '0', color='black', fontsize=4)
#n3graph.grid(True)
#The following code subplots the I1/I3 vs time graph
i1i3ticks = np.arange(-5,75,5)
i1i3graph.plot(data['time'], data['freqI1I3'])
#i1i3graph.set_title('I1/I3 Input', fontsize = 4)
i1i3graph.set_xlabel('Time (seconds)', fontsize = 2)
i1i3graph.set_ylabel('I1/I3 Frequency (Hz)', fontsize = 2)
i1i3graph.axis([0, 9, -5, 75])
i1i3graph.set_xticks(timeticks)
i1i3graph.set_yticks(i1i3ticks)
i1i3graph.tick_params(labelsize=3)
i1i3graph.set_axis_off()
i1i3graph.text(9.0140,10, '- I 1 I 3 -', color='black', fontsize=4)
i1i3graph.text(-0.5140,70, '70', color='black', fontsize=4)
i1i3graph.text(-0.5140,00, '0', color='black', fontsize=4)
#i1i3graph.grid(True)

#The following code subplots the Hinge input vs time graph
hingeticks = np.arange(-5,25, 5)
hingegraph.plot(data['time'], data['freqHinge'])
#hingegraph.set_title('Hinge Input', fontsize = 4)
#hingegraph.set_xlabel('Time (seconds)', fontsize = 2)
hingegraph.set_ylabel('Hinge Frequency (Hz)', fontsize = 2)
hingegraph.axis([0, 9, -5, 25], fontsize = 2)
hingegraph.set_xticks(timeticks)
hingegraph.set_yticks(hingeticks)
hingegraph.tick_params(labelsize=3)
hingegraph.set_axis_off()
hingegraph.text(9.0140,10, '- HINGE -', color='black', fontsize=4)
hingegraph.text(-0.5140,20, '20', color='black', fontsize=4)
hingegraph.text(-0.5140,00, '0', color='black', fontsize=4)
#hingegraph.grid(True)

vertlinei2 = i2graph.vlines([data['time'][0]],[0],80, colors = 'r')
vertlinen3 = n3graph.vlines([data['time'][0]],[0],80, colors = 'r')
vertlinei1i3 = i1i3graph.vlines([data['time'][0]],[0],80, colors = 'r')
vertlinehinge = hingegraph.vlines([data['time'][0]],[0],80, colors = 'r')

def makevertlinei2(i2graph, i):
    return i2graph.vlines([data['time'][i]],[0],80, colors = 'r')

def makevertlinen3(n3graph, i):
    return n3graph.vlines([data['time'][i]],[0],80, colors = 'r')

def makevertlinei1i3(i1i3graph, i):
    return i1i3graph.vlines([data['time'][i]],[0],80, colors = 'r')

def makevertlinehinge(hingegraph, i):
    return hingegraph.vlines([data['time'][i]],[0],80, colors = 'r')

def createOdontophore(i, odontophore):
    #Define Variables by acsessing array values
    oradius = data['a'][i] #"a" in the kinetic model
    odiameter = 2*oradius #2*a
    ominoraxisradius = (5.0 * sqrt(5.0) / sqrt(oradius * 1000.0)) / 1000.0 #"b" in the kinetic model
    ominoraxisdiameter = 2*ominoraxisradius #2*b
    oangle = -(90 - data['odontophoreangle'][i]) #"odontophoreangle"
    ocenter = data['x'][i] #"x"
    #odontophore
    odontophore = Ellipse((ocenter, 0), ominoraxisdiameter, odiameter, oangle)
    odontophore.set_alpha(1)
    odontophore.set_facecolor('none')
    odontophore.set_linewidth(.5)
    if data['freqN3'][i] == 0:
        odontophore.set_edgecolor('black')
    elif data['freqN3'][i] > 0:
        odontophore.set_edgecolor('black')
    else:
        odontophore.set_edgecolor('red')
    return odontophore

def createi1i3top(i, i1i3top):
    topcontactpointy = data['ytop'][i]
    i1i3top = Circle((0,topcontactpointy),0.00125)
    i1i3top.set_facecolor('none')
    i1i3top.set_linewidth(.5)
    if data['freqI1I3'][i] == 0:
        i1i3top.set_edgecolor('black')
    elif data['freqI1I3'][i] > 0:
        i1i3top.set_edgecolor('black')
    else:
        i1i3top.set_edgecolor('red')
    return i1i3top

def createi1i3bottom(i, i1i3bottom):
    bottomcontactpointy = data['ybottom'][i]
    i1i3bottom = Circle((0,bottomcontactpointy),0.00125)
    i1i3bottom.set_facecolor('none')
    i1i3bottom.set_linewidth(.5)
    if data['freqI1I3'][i] == 0:
        i1i3bottom.set_edgecolor('black')
    elif data['freqI1I3'][i] > 0:
        i1i3bottom.set_edgecolor('black')
    else:
        i1i3bottom.set_edgecolor('red')
    return i1i3bottom

def createi2top(i, i2top):
    #Top line

    i2length = data['lengthofI2'][i]
    topangle = data['topphiangleofi2'][i] #Tphi
    bottomangle = data['bottomphiangleofi2'][i] #Bphi
    xtopval = -(np.cos(np.deg2rad(topangle))*i2length)/2
    xbottomval = -(np.cos(np.deg2rad(bottomangle))*i2length)/2

    furthestbackxpoint = data['furthestbackxpoint'][i]
    furthestbackypoint = data['furthestbackypoint'][i]
    i1i3contacttopy = data['i1i3contacttopy'][i]
    i1i3contacttopx = -0.00125
    i1i3contactbottomy = data['i1i3contactbottomy'][i]
    i1i3contactbottomx = -.00125

    ocontacttopx = data['ocontacttopx'][i]
    ocontacttopy = data['ocontacttopy'][i]
    ocontactbottomx = data['ocontactbottomx'][i]
    ocontactbottomy = data['ocontactbottomy'][i]

    x1,y1 = np.array([[i1i3contacttopx, ocontacttopx],[i1i3contacttopy, ocontacttopy]])
    i2top = mlines.Line2D(x1,y1)
    i2top.set_linewidth(.5)
    if data['freqI2'][i] == 0:
        i2top.set_color('black')
    elif data['freqI2'][i] > 0:
        i2top.set_color('black')
    else:
        i2top.set_color('red')
    return i2top

def createi2bottom(i, i2bottom):
    #Bottom line

    i2length = data['lengthofI2'][i]
    topangle = data['topphiangleofi2'][i] #Tphi
    bottomangle = data['bottomphiangleofi2'][i] #Bphi
    xtopval = -(np.cos(np.deg2rad(topangle))*i2length)/2
    xbottomval = -(np.cos(np.deg2rad(bottomangle))*i2length)/2

    furthestbackxpoint = data['furthestbackxpoint'][i]
    furthestbackypoint = data['furthestbackypoint'][i]
    i1i3contacttopy = data['i1i3contacttopy'][i]
    i1i3contacttopx = -0.00125
    i1i3contactbottomy = data['i1i3contactbottomy'][i]
    i1i3contactbottomx = -.00125

    ocontacttopx = data['ocontacttopx'][i]
    ocontacttopy = data['ocontacttopy'][i]
    ocontactbottomx = data['ocontactbottomx'][i]
    ocontactbottomy = data['ocontactbottomy'][i]

    x2, y2 = np.array([[i1i3contactbottomx, ocontactbottomx],[i1i3contactbottomy, ocontactbottomy]])
    i2bottom= mlines.Line2D(x2,y2)
    i2bottom.set_linewidth(.5)
    if data['freqI2'][i] == 0:
        i2bottom.set_color('black')
    elif data['freqI2'][i] > 0:
        i2bottom.set_color('black')
    else:
        i2bottom.set_color('red')
    return i2bottom

def maketoptextI2(figure, i, toptexti2):
    if data['freqI2'][i] == 0:
        toptexti2 = figure.text(.0140,.0175, '- I 2 -', color='black', fontsize=4)
    elif data['freqI2'][i] > 0:
        toptexti2 = figure.text(.0140,.0175, '- I 2 -', color='green', fontsize=4)
    else:
        toptexti2 = figure.text(.0140,.0175, '%f' % data['freqI2'][i], color='red', fontsize=4)
    return toptexti2

def maketoptextN3(figure, i, toptextN3):
    if data['freqN3'][i] == 0:
        toptextN3 = figure.text(-.008,.0175, '- I 4 -', color='black', fontsize=4)
    elif data['freqN3'][i] > 0:
        toptextN3 = figure.text(-.008,.0175, '- I 4 -', color='green', fontsize=4)
    else:
        toptextN3 = figure.text(-.008,.0175, '%f' % data['freqN3'][i], color='red', fontsize=4)
    return toptextN3

def maketoptexti1i3(figure, i, toptexti1i3):
    if data['freqI1I3'][i] == 0:
        toptexti1i3 = figure.text(-.0020,.0175, '- I 1 I 3 -', color='black', fontsize=4)
    elif data['freqI1I3'][i] > 0:
        toptexti1i3 = figure.text(-.0020,.0175, '- I 1 I 3 -', color='green', fontsize=4)
    else:
        toptexti1i3 = figure.text(-.0020,.0175, '%f' % data['freqI1I3'][i], color='red', fontsize=4)
    return toptexti1i3

def maketoptexthinge(figure, i, toptexthinge):
    if data['freqHinge'][i] == 0:
        toptexthinge = figure.text(.0045,.0175, '- H I N G E -', color='black', fontsize=4)
    elif data['freqHinge'][i] > 0:
        toptexthinge = figure.text(.0045,.0175, '- H I N G E -', color='green', fontsize=4)
    else:
        toptexthinge = figure.text(.0025,.0175, '%f' % data['freqHinge'][i], color='red', fontsize=4)
    return toptexthinge

def makebottomtexti2(figure, i, bottomtexti2):
    if data['freqI2'][i] == 0:
        bottomtexti2 = figure.text(.0140,.015, '%2.0f' % data['freqI2'][i], color='black', fontsize=4)
    elif data['freqI2'][i] > 0:
        bottomtexti2 = figure.text(.0140,.015, '%2.0f' % data['freqI2'][i], color='green', fontsize=4)
    else:
        bottomtexti2 = figure.text(.0140,.015, '%f' % data['freqI2'][i], color='red', fontsize=4)
    return bottomtexti2

def makebottomtextN3(figure, i, bottomtextN3):
    if data['freqN3'][i] == 0:
        bottomtextN3 = figure.text(-.008,.015, '%2.0f' % data['freqN3'][i], color='black', fontsize=4)
    elif data['freqN3'][i] > 0:
        bottomtextN3 = figure.text(-.008,.015, '%2.0f' % data['freqN3'][i], color='green', fontsize=4)
    else:
        bottomtextN3 = figure.text(-.008,.015, '%f' % data['freqN3'][i], color='red', fontsize=4)
    return bottomtextN3

def makebottomtexti1i3(figure, i, bottomtexti1i3):
    if data['freqI1I3'][i] == 0:
        bottomtexti1i3 = figure.text(-.0020,.015, '%2.0f' % data['freqI1I3'][i], color='black', fontsize=4)
    elif data['freqI1I3'][i] > 0:
        bottomtexti1i3 = figure.text(-.0020,.015, '%2.0f' % data['freqI1I3'][i], color='green', fontsize=4)
    else:
        bottomtexti1i3 = figure.text(-.0020,.015, '%f' % data['freqI1I3'][i], color='red', fontsize=4)
    return bottomtexti1i3

def makebottomtexthinge(figure, i, bottomtexthinge):
    if data['freqHinge'][i] == 0:
        bottomtexthinge = figure.text(.0045,.015, '%2.0f' % data['freqHinge'][i], color='black', fontsize=4)
    elif data['freqHinge'][i] > 0:
        bottomtexthinge = figure.text(.0045,.015, '%2.0f' % data['freqHinge'][i], color='green', fontsize=4)
    else:
        bottomtexthinge = figure.text(.0045,.015, '%f' % data['freqHinge'][i], color='red', fontsize=4)
    return bottomtexthinge

def maketimer(figure, i, timer):
    timer = figure.text(0, .02, 'Time: %3.2f seconds' % float(i/100.00), color='Black', fontsize=6, horizontalalignment = 'center')
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
elif animation.AVConvWriter.isAvailable():
    movieWriterClass = animation.AVConvWriter
else:
    sys.stderr.write('video converter missing')
    exit()

#For each i/time step, save the figure as a frame in the .mp4 file
moviewriter = movieWriterClass(fps = fps, bitrate = br)
with moviewriter.saving(plt.gcf(), moviefile, dpi):
    newFrameEveryNLines = 5  # 5 * 0.01 sec = 0.05 sec sim time per frame = 1/fps = playback in realtime
    for i in tqdm(np.arange(0, data.shape[0], newFrameEveryNLines)):
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
