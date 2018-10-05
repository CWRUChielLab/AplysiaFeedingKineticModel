import os
import sys
import csv
import pandas as pd
import matplotlib
matplotlib.use('Agg') # force matplotlib to not use XWindows backend
from math import sqrt
import matplotlib.gridspec as gridspec
from scipy.interpolate import spline
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.lines as mlines
from matplotlib.patches import Circle
from matplotlib import animation

class Animator:
    def __init__(self):
        self.setVideoParameters()
        self.readFile()
        self.createAnimationBox()
        self.initShapeObjects()
        self.createMovie()

    def saveMovie(self):
        self.createMovie()
        
######################################################
##                 Video Parameters                 ##
######################################################
    def setVideoParameters(self):
        #File Name
        self.moviefile = "animation.mp4"
        #size of movie
        self.movieWidth = 1600
        self.height = 800
        #frames per second
        self.fps = 20
        #bitrate to reduce compression
        self.br = 2000
        #moviesize
        self.dpi = 180

######################################################
##                    Read File                     ##
######################################################
    def readFile(self):
        self.data = pd.read_csv('output-mechanics.csv') 

######################################################
##                Plot 2D Simulation                ##
######################################################
    def createAnimationBox(self):
        self.figure = plt.figure()
        self.gs = gridspec.GridSpec(18, 9)
        self.a = self.figure.add_subplot(self.gs[2:9,0:9], aspect='equal')
        plt.ylabel('Y- Plane (m)', fontsize = 4)
        plt.xlabel('X- Plane (m)', fontsize = 4)
        plt.xlim(-.02,.02)
        plt.ylim(-.01,.01)
        plt.title('Animation', fontsize = 8)


    ##Jeff, this is currently not used but can be added
    def createFreqPlots(self):
        # graphs for frequencies, a is animation figure
        self.i2graph = self.figure.add_subplot(gs[11:12,1:8])
        self.n3graph = self.figure.add_subplot(gs[13:14,1:8])
        self.i1i3graph = self.figure.add_subplot(gs[15:16,1:8])
        self.hingegraph = self.figure.add_subplot(gs[17:18,1:8])

    def initShapeObjects(self):
        # init shapes
        self.odontophore = Ellipse((0, 0), 0, 0, 0)
        self.i1i3top = Circle((0,0),0.00125)
        self.i1i3bottom = Circle((0,0),0.00125)
        self.i1i3top1 = Circle((0,0),0.00125)
        self.i1i3bottom1 = Circle((0,0),0.00125)
        self.i1i3top2 = Circle((0,0),0.00125)
        self.i1i3bottom2 = Circle((0,0),0.00125)
        self.i1i3top3 = Circle((0,0),0.00125)
        self.i1i3bottom3 = Circle((0,0),0.00125)
        self.i1i3top4 = Circle((0,0),0.00125)
        self.i1i3bottom4 = Circle((0,0),0.00125)
        self.i2top = mlines.Line2D([0,1],[0,1])
        self.i2bottom= mlines.Line2D([0,1],[0,1])
        # add them to the animation box
        self.a.add_artist(self.odontophore)
        self.a.add_artist(self.i1i3top)
        self.a.add_artist(self.i1i3bottom)
        self.a.add_artist(self.i1i3top1)
        self.a.add_artist(self.i1i3bottom1)
        self.a.add_artist(self.i1i3top2)
        self.a.add_artist(self.i1i3bottom2)
        self.a.add_artist(self.i1i3top3)
        self.a.add_artist(self.i1i3bottom3)
        self.a.add_artist(self.i1i3top4)
        self.a.add_artist(self.i1i3bottom4)
        self.a.add_artist(self.i2top)
        self.a.add_artist(self.i2bottom)
        #plt.ylabel('Y- Plane (m)', fontsize = 4)
        #plt.xlabel('X- Plane (m)', fontsize = 4)
        #plt.xlim(-.02,.02)
        #plt.ylim(-.01,.01)
        #plt.title('Animation', fontsize = 8)

    ##Jeff, this is currently not used but can be added
    def createFreqPlots(self):
        self.a.tick_params(labelsize=3)
        self.a.text(-.0240,.0175, 'Channels Activated:', color='Black', fontsize=4)
        self.a.text(-.0240,.015, 'Level of Activation (Hz):', color='Black', fontsize=4)
        #a text
        self.toptextI2 = a.text(.0140,.0175, '- I 2 -', color='black', fontsize=4)
        self.toptextN3 = a.text(-.008,.0175, '- I 4 -', color='black', fontsize=4)
        self.toptexti1i3 = a.text(-.0020,.0175, '- I 1 I 3 -', color='black', fontsize=4)
        self.toptexthinge = a.text(.0045,.0175, '- H I N G E -', color='black', fontsize=4)
        self.bottomtexti2 = a.text(.0140,.015, '%2.1f' % float(freqi2array[1]), color='black', fontsize=4)
        self.bottomtextN3 = a.text(-.008,.015, '%2.1f' % float(freqN3array[1]), color='black', fontsize=4)
        self.bottomtexti1i3 = a.text(-.0020,.015, '%2.1f' % float(freqi1i3array[1]), color='black', fontsize=4)
        self.bottomtexthinge = a.text(.0045,.015, '%2.1f' % float(freqHingearray[1]), color='black', fontsize=4)
        self.timer = a.text(0, .02, 'Time: 0 seconds', color='Black', fontsize=6, horizontalalignment = 'center')

        self.timeticks = np.arange(0,9,1)
        #The following code subplots the I2 input vs time graph
        self.i2ticks = np.arange(-5,35,5)
        self.i2graph.plot(time[1:850],freqi2array[1:850])
        #i2graph.set_title('I2 Input', fontsize = 4)
        #i2graph.set_xlabel(time[0] + ' (seconds)', fontsize = 2)
        self.i2graph.set_ylabel(freqi2array[0] + ' (Hz)', fontsize = 2)
        self.i2graph.axis([0, 9, -5, 35])
        self.i2graph.set_xticks(timeticks)
        self.i2graph.set_yticks(i2ticks)
        self.i2graph.tick_params(labelsize=3)
        self.i2graph.set_axis_off()
        self.i2graph.text(9.0140,10, '- I 2 -', color='black', fontsize=4)
        self.i2graph.text(-0.5140,30, '30', color='black', fontsize=4)
        self.i2graph.text(-0.5140,00, '0', color='black', fontsize=4)
        #i2graph.grid(True)
            
        #The following code subplots the Odontophore input vs time graph
        self.n3ticks = np.arange(-5,40, 5)
        self.n3graph.plot(time[1:850],freqN3array[1:850])
        #n3graph.set_title('Odontophore Input', fontsize = 4)
        #n3graph.set_xlabel(time[0] + ' (seconds)', fontsize = 2)
        self.n3graph.set_ylabel(freqN3array[0] + ' (Hz)', fontsize = 2)
        self.n3graph.axis([0, 9, -5, 40], fontsize = 2)
        self.n3graph.set_xticks(timeticks)
        self.n3graph.set_yticks(n3ticks)
        self.n3graph.tick_params(labelsize=3)
        self.n3graph.set_axis_off()
        self.n3graph.text(9.0140,10, '- I 4 -', color='black', fontsize=4)
        self.n3graph.text(-0.5140,35, '35', color='black', fontsize=4)
        self.n3graph.text(-0.5140,00, '0', color='black', fontsize=4)
        #n3graph.grid(True)
        #The following code subplots the I1/I3 vs time graph
        self.i1i3ticks = np.arange(-5,75,5)
        self.i1i3graph.plot(time[1:850],freqi1i3array[1:850])
        #i1i3graph.set_title('I1/I3 Input', fontsize = 4)
        self.i1i3graph.set_xlabel(time[0] + ' (seconds)', fontsize = 2)
        self.i1i3graph.set_ylabel(freqi1i3array[0] + ' (Hz)', fontsize = 2)
        self.i1i3graph.axis([0, 9, -5, 75])
        self.i1i3graph.set_xticks(timeticks)
        self.i1i3graph.set_yticks(i1i3ticks)
        self.i1i3graph.tick_params(labelsize=3)
        self.i1i3graph.set_axis_off()
        self.i1i3graph.text(9.0140,10, '- I 1 I 3 -', color='black', fontsize=4)
        self.i1i3graph.text(-0.5140,70, '70', color='black', fontsize=4)
        self.i1i3graph.text(-0.5140,00, '0', color='black', fontsize=4)
        #i1i3graph.grid(True)

        #The following code subplots the Hinge input vs time graph
        self.hingeticks = np.arange(-5,25, 5)
        self.hingegraph.plot(time[1:850],freqHingearray[1:850])
        #hingegraph.set_title('Hinge Input', fontsize = 4)
        #hingegraph.set_xlabel(time[0] + ' (seconds)', fontsize = 2)
        self.hingegraph.set_ylabel(freqHingearray[0] + ' (Hz)', fontsize = 2)
        self.hingegraph.axis([0, 9, -5, 25], fontsize = 2)
        self.hingegraph.set_xticks(timeticks)
        self.hingegraph.set_yticks(hingeticks)
        self.hingegraph.tick_params(labelsize=3)
        self.hingegraph.set_axis_off()
        self.hingegraph.text(9.0140,10, '- HINGE -', color='black', fontsize=4)
        self.hingegraph.text(-0.5140,20, '20', color='black', fontsize=4)
        self.hingegraph.text(-0.5140,00, '0', color='black', fontsize=4)
        #hingegraph.grid(True)

        self.vertlinei2 = i2graph.vlines([float(time[1])],[0],80, colors = 'r')
        self.vertlinen3 = n3graph.vlines([float(time[1])],[0],80, colors = 'r')
        self.vertlinei1i3 = i1i3graph.vlines([float(time[1])],[0],80, colors = 'r')
        self.vertlinehinge = hingegraph.vlines([float(time[1])],[0],80, colors = 'r')

    def makevertlinei2(self, i2graph, discretemoment):
        return self.i2graph.vlines([float(time[discretemoment])],[0],80, colors = 'r')
        
    def makevertlinen3(self, n3graph, discretemoment):
        return self.n3graph.vlines([float(time[discretemoment])],[0],80, colors = 'r')
        
    def makevertlinei1i3(self, i1i3graph, discretemoment):
        return self.i1i3graph.vlines([float(time[discretemoment])],[0],80, colors = 'r')
        
    def makevertlinehinge(self, hingegraph, discretemoment):
        return self.hingegraph.vlines([float(time[discretemoment])],[0],80, colors = 'r')

    def createOdontophore(self, discretemoment, odontophore):
        #Define Variables by acsessing array values
        self.oradius = float(self.data['a'][discretemoment]) #"a" in the kinetic model
        self.odiameter = 2*self.oradius #2*a
        self.ominoraxisradius = (5.0 * sqrt(5.0) / sqrt(self.oradius * 1000.0)) / 1000.0 #"b" in the kinetic model
        self.ominoraxisdiameter = 2*self.ominoraxisradius #2*b
        self.oangle = -(90 - float(self.data['odontophoreangle'][discretemoment])) #"odontophoreangle"
        self.ocenter = float(self.data['x'][discretemoment]) #"x"
        #odontophore
        self.odontophore = Ellipse((self.ocenter, 0), self.ominoraxisdiameter, self.odiameter, self.oangle)
        self.odontophore.set_alpha(1)
        self.odontophore.set_facecolor('none')
        self.odontophore.set_linewidth(.5)
        self.odontophore.set_edgecolor('black')
    ##    if float(freqN3array[discretemoment]) == 0:
    ##        self.odontophore.set_edgecolor('black')
    ##    elif float(freqN3array[discretemoment]) > 0:
    ##        self.odontophore.set_edgecolor('black')
    ##    else:
    ##        self.odontophore.set_edgecolor('red')
        return self.odontophore

    def createi1i3top(self, discretemoment, i1i3top):
        self.topcontactpointy = float(self.data['ytop'][discretemoment])
        self.i1i3top = Circle((0,self.topcontactpointy),0.00125)
        self.i1i3top.set_facecolor('none')
        self.i1i3top.set_linewidth(.5)
        self.i1i3top.set_edgecolor('black')

    ##    if float(freqi1i3array[discretemoment]) == 0:
    ##        self.i1i3top.set_edgecolor('black')
    ##    elif float(freqi1i3array[discretemoment]) > 0:
    ##        self.i1i3top.set_edgecolor('black')
    ##    else:
    ##        self.i1i3top.set_edgecolor('red')
        return self.i1i3top

    def createi1i3bottom(self, discretemoment, i1i3bottom):
        self.bottomcontactpointy = float(self.data['ybottom'][discretemoment])
        self.i1i3bottom = Circle((0,self.bottomcontactpointy),0.00125)
        self.i1i3bottom.set_facecolor('none')
        self.i1i3bottom.set_linewidth(.5)
        self.i1i3bottom.set_edgecolor('black')
    ##    if float(freqi1i3array[discretemoment]) == 0:
    ##        self.i1i3bottom.set_edgecolor('black')
    ##    elif float(freqi1i3array[discretemoment]) > 0:
    ##        self.i1i3bottom.set_edgecolor('black')
    ##    else:
    ##        self.i1i3bottom.set_edgecolor('red')
        return self.i1i3bottom

    def createi1i3top1(self, discretemoment, i1i3top1):
        self.topcontactpointy1 = float(self.data['ytop1'][discretemoment])
        self.i1i3top1 = Circle((1*2*0.00125,self.topcontactpointy1),0.00125)
        self.i1i3top1.set_facecolor('none')
        self.i1i3top1.set_linewidth(.5)
        self.i1i3top1.set_edgecolor('black')

    ##    if float(freqi1i3array[discretemoment]) == 0:
    ##        self.i1i3top.set_edgecolor('black')
    ##    elif float(freqi1i3array[discretemoment]) > 0:
    ##        self.i1i3top.set_edgecolor('black')
    ##    else:
    ##        self.i1i3top.set_edgecolor('red')
        return self.i1i3top1

    def createi1i3bottom1(self, discretemoment, i1i3bottom1):
        self.bottomcontactpointy1 = float(self.data['ybottom1'][discretemoment])
        self.i1i3bottom1 = Circle((1*2*0.00125,self.bottomcontactpointy1),0.00125)
        self.i1i3bottom1.set_facecolor('none')
        self.i1i3bottom1.set_linewidth(.5)
        self.i1i3bottom1.set_edgecolor('black')
    ##    if float(freqi1i3array[discretemoment]) == 0:
    ##        self.i1i3bottom.set_edgecolor('black')
    ##    elif float(freqi1i3array[discretemoment]) > 0:
    ##        self.i1i3bottom.set_edgecolor('black')
    ##    else:
    ##        self.i1i3bottom.set_edgecolor('red')
        return self.i1i3bottom1

    def createi1i3top2(self, discretemoment, i1i3top2):
        self.topcontactpointy2 = float(self.data['ytop2'][discretemoment])
        self.i1i3top2 = Circle((2*2*0.00125,self.topcontactpointy2),0.00125)
        self.i1i3top2.set_facecolor('none')
        self.i1i3top2.set_linewidth(.5)
        self.i1i3top2.set_edgecolor('black')

    ##    if float(freqi1i3array[discretemoment]) == 0:
    ##        self.i1i3top.set_edgecolor('black')
    ##    elif float(freqi1i3array[discretemoment]) > 0:
    ##        self.i1i3top.set_edgecolor('black')
    ##    else:
    ##        self.i1i3top.set_edgecolor('red')
        return self.i1i3top2

    def createi1i3bottom2(self, discretemoment, i1i3bottom2):
        self.bottomcontactpointy2 = float(self.data['ybottom1'][discretemoment])
        self.i1i3bottom2 = Circle((2*2*0.00125,self.bottomcontactpointy2),0.00125)
        self.i1i3bottom2.set_facecolor('none')
        self.i1i3bottom2.set_linewidth(.5)
        self.i1i3bottom2.set_edgecolor('black')
    ##    if float(freqi1i3array[discretemoment]) == 0:
    ##        self.i1i3bottom.set_edgecolor('black')
    ##    elif float(freqi1i3array[discretemoment]) > 0:
    ##        self.i1i3bottom.set_edgecolor('black')
    ##    else:
    ##        self.i1i3bottom.set_edgecolor('red')
        return self.i1i3bottom2

    def createi1i3top3(self, discretemoment, i1i3top3):
        self.topcontactpointy3 = float(self.data['ytop3'][discretemoment])
        self.i1i3top3 = Circle((3*2*0.00125,self.topcontactpointy3),0.00125)
        self.i1i3top3.set_facecolor('none')
        self.i1i3top3.set_linewidth(.5)
        self.i1i3top3.set_edgecolor('black')

    ##    if float(freqi1i3array[discretemoment]) == 0:
    ##        self.i1i3top.set_edgecolor('black')
    ##    elif float(freqi1i3array[discretemoment]) > 0:
    ##        self.i1i3top.set_edgecolor('black')
    ##    else:
    ##        self.i1i3top.set_edgecolor('red')
        return self.i1i3top3

    def createi1i3bottom3(self, discretemoment, i1i3bottom3):
        self.bottomcontactpointy3 = float(self.data['ybottom3'][discretemoment])
        self.i1i3bottom3 = Circle((3*2*0.00125,self.bottomcontactpointy3),0.00125)
        self.i1i3bottom3.set_facecolor('none')
        self.i1i3bottom3.set_linewidth(.5)
        self.i1i3bottom3.set_edgecolor('black')
    ##    if float(freqi1i3array[discretemoment]) == 0:
    ##        self.i1i3bottom.set_edgecolor('black')
    ##    elif float(freqi1i3array[discretemoment]) > 0:
    ##        self.i1i3bottom.set_edgecolor('black')
    ##    else:
    ##        self.i1i3bottom.set_edgecolor('red')
        return self.i1i3bottom3

    def createi1i3top4(self, discretemoment, i1i3top4):
        self.topcontactpointy4 = float(self.data['ytop4'][discretemoment])
        self.i1i3top4 = Circle((4*2*0.00125,self.topcontactpointy4),0.00125)
        self.i1i3top4.set_facecolor('none')
        self.i1i3top4.set_linewidth(.5)
        self.i1i3top4.set_edgecolor('black')

    ##    if float(freqi1i3array[discretemoment]) == 0:
    ##        self.i1i3top.set_edgecolor('black')
    ##    elif float(freqi1i3array[discretemoment]) > 0:
    ##        self.i1i3top.set_edgecolor('black')
    ##    else:
    ##        self.i1i3top.set_edgecolor('red')
        return self.i1i3top4

    def createi1i3bottom4(self, discretemoment, i1i3bottom4):
        self.bottomcontactpointy4 = float(self.data['ybottom4'][discretemoment])
        self.i1i3bottom4 = Circle((4*2*0.00125,self.bottomcontactpointy4),0.00125)
        self.i1i3bottom4.set_facecolor('none')
        self.i1i3bottom4.set_linewidth(.5)
        self.i1i3bottom4.set_edgecolor('black')
    ##    if float(freqi1i3array[discretemoment]) == 0:
    ##        self.i1i3bottom.set_edgecolor('black')
    ##    elif float(freqi1i3array[discretemoment]) > 0:
    ##        self.i1i3bottom.set_edgecolor('black')
    ##    else:
    ##        self.i1i3bottom.set_edgecolor('red')
        return self.i1i3bottom4


    def createi2top(self, discretemoment, i2top):
        #Top line
        
        #self.i1i3radius = float(i1i3radiusarray[discretemoment])
        self.i2length = float(self.data['lengthofI2'][discretemoment])
        self.topangle = float(self.data['topphiangleofi2'][discretemoment]) #Tphi
        self.bottomangle = float(self.data['bottomphiangleofi2'][discretemoment]) #Bphi
        self.xtopval = -(np.cos(np.deg2rad(self.topangle))*self.i2length)/2
        self.xbottomval = -(np.cos(np.deg2rad(self.bottomangle))*self.i2length)/2
        
        #self.furthestbackxpoint = float(furthestbackxpointarray[discretemoment])
        #self.furthestbackypoint = float(furthestbackypointarray[discretemoment])
        self.i1i3contacttopy = float(self.data['i1i3contacttopy'][discretemoment])
        self.i1i3contacttopx = -0.00125
        self.i1i3contactbottomy = float(self.data['i1i3contactbottomy'][discretemoment])
        self.i1i3contactbottomx = -.00125
        
        #x1val = float(x1array[discretemoment])
        
        self.ocontacttopx = float(self.data['ocontacttopx'][discretemoment])
        self.ocontacttopy = float(self.data['ocontacttopy'][discretemoment])
        self.ocontactbottomx = float(self.data['ocontactbottomx'][discretemoment])
        self.ocontactbottomy = float(self.data['ocontactbottomy'][discretemoment])

        x1,y1 = np.array([[self.i1i3contacttopx, self.ocontacttopx],[self.i1i3contacttopy, self.ocontacttopy]])
        self.i2top = mlines.Line2D(x1,y1)
        self.i2top.set_linewidth(.5)
        self.i2top.set_color('black')
    ##    if float(freqi2array[discretemoment]) == 0:
    ##        self.i2top.set_color('black')
    ##    elif float(freqi2array[discretemoment]) > 0:
    ##        self.i2top.set_color('black')
    ##    else:
    ##        self.i2top.set_color('red')
        return self.i2top

    def createi2bottom(self, discretemoment, i2bottom):
        #Bottom line
        
        #self.i1i3radius = float(i1i3radiusarray[discretemoment])
        self.i2length = float(self.data['lengthofI2'][discretemoment])
        self.topangle = float(self.data['topphiangleofi2'][discretemoment]) #Tphi
        self.bottomangle = float(self.data['bottomphiangleofi2'][discretemoment]) #Bphi
        self.xtopval = -(np.cos(np.deg2rad(self.topangle))*self.i2length)/2
        self.xbottomval = -(np.cos(np.deg2rad(self.bottomangle))*self.i2length)/2
        
        #self.furthestbackxpoint = float(furthestbackxpointarray[discretemoment])
        #self.furthestbackypoint = float(furthestbackypointarray[discretemoment])
        self.i1i3contacttopy = float(self.data['i1i3contacttopy'][discretemoment])
        self.i1i3contacttopx = -0.00125
        self.i1i3contactbottomy = float(self.data['i1i3contactbottomy'][discretemoment])
        self.i1i3contactbottomx = -.00125
        
        #x1val = float(x1array[discretemoment])
        
        self.ocontacttopx = float(self.data['ocontacttopx'][discretemoment])
        self.ocontacttopy = float(self.data['ocontacttopy'][discretemoment])
        self.ocontactbottomx = float(self.data['ocontactbottomx'][discretemoment])
        self.ocontactbottomy = float(self.data['ocontactbottomy'][discretemoment])

        x2, y2 = np.array([[self.i1i3contactbottomx, self.ocontactbottomx],[self.i1i3contactbottomy,self.ocontactbottomy]])
        self.i2bottom= mlines.Line2D(x2,y2)
        self.i2bottom.set_linewidth(.5)
        self.i2bottom.set_color('black')
    ##    if float(freqi2array[discretemoment]) == 0:
    ##        self.i2bottom.set_color('black')
    ##    elif float(freqi2array[discretemoment]) > 0:
    ##        self.i2bottom.set_color('black')
    ##    else:
    ##        self.i2bottom.set_color('red')
        return self.i2bottom

    def maketoptextI2(self, figure, discretemoment, toptexti2):
        if float(freqi2array[discretemoment]) == 0:
            self.toptexti2 = figure.text(.0140,.0175, '- I 2 -', color='black', fontsize=4)
        elif float(freqi2array[discretemoment]) > 0:
            self.toptexti2 = figure.text(.0140,.0175, '- I 2 -', color='green', fontsize=4)
        else:
            self.toptexti2 = figure.text(.0140,.0175, '%f' % float(freqi2array[discretemoment]), color='red', fontsize=4)
        return self.toptexti2

    def maketoptextN3(figure, discretemoment, toptextN3):
        if float(freqN3array[discretemoment]) == 0:
            self.toptextN3 = figure.text(-.008,.0175, '- I 4 -', color='black', fontsize=4)
        elif float(freqN3array[discretemoment]) > 0:
            self.toptextN3 = figure.text(-.008,.0175, '- I 4 -', color='green', fontsize=4)
        else:
            self.toptextN3 = figure.text(-.008,.0175, '%f' % float(freqN3array[discretemoment]), color='red', fontsize=4)
        return self.toptextN3

    def maketoptexti1i3(figure, discretemoment, toptexti1i3):
        if float(freqi1i3array[discretemoment]) == 0:
            self.toptexti1i3 = figure.text(-.0020,.0175, '- I 1 I 3 -', color='black', fontsize=4)
        elif float(freqi1i3array[discretemoment]) > 0:
            self.toptexti1i3 = figure.text(-.0020,.0175, '- I 1 I 3 -', color='green', fontsize=4)
        else:
            self.toptexti1i3 = figure.text(-.0020,.0175, '%f' % float(freqi1i3array[discretemoment]), color='red', fontsize=4)
        return self.toptexti1i3

    def maketoptexthinge(figure, discretemoment, toptexthinge):
        if float(freqHingearray[discretemoment]) == 0:
            self.toptexthinge = figure.text(.0045,.0175, '- H I N G E -', color='black', fontsize=4)
        elif float(freqHingearray[discretemoment]) > 0:
            self.toptexthinge = figure.text(.0045,.0175, '- H I N G E -', color='green', fontsize=4)
        else:
            self.toptexthinge = figure.text(.0025,.0175, '%f' % float(freqHingearray[discretemoment]), color='red', fontsize=4)
        return self.toptexthinge

    def makebottomtexti2(figure, discretemoment, bottomtexti2):
        if float(freqi2array[discretemoment]) == 0:
            self.bottomtexti2 = figure.text(.0140,.015, '%2.0f' % float(freqi2array[discretemoment]), color='black', fontsize=4)
        elif float(freqi2array[discretemoment]) > 0:
            self.bottomtexti2 = figure.text(.0140,.015, '%2.0f' % float(freqi2array[discretemoment]), color='green', fontsize=4)
        else:
            self.bottomtexti2 = figure.text(.0140,.015, '%f' % float(freqi2array[discretemoment]), color='red', fontsize=4)
        return self.bottomtexti2

    def makebottomtextN3(figure, discretemoment, bottomtextN3):
        if float(freqN3array[discretemoment]) == 0:
            self.bottomtextN3 = figure.text(-.008,.015, '%2.0f' % float(freqN3array[discretemoment]), color='black', fontsize=4)
        elif float(freqN3array[discretemoment]) > 0:
            self.bottomtextN3 = figure.text(-.008,.015, '%2.0f' % float(freqN3array[discretemoment]), color='green', fontsize=4)
        else:
            self.bottomtextN3 = figure.text(-.008,.015, '%f' % float(freqN3array[discretemoment]), color='red', fontsize=4)
        return self.bottomtextN3

    def makebottomtexti1i3(figure, discretemoment, bottomtexti1i3):
        if float(freqi1i3array[discretemoment]) == 0:
            self.bottomtexti1i3 = figure.text(-.0020,.015, '%2.0f' % float(freqi1i3array[discretemoment]), color='black', fontsize=4)
        elif float(freqi1i3array[discretemoment]) > 0:
            self.bottomtexti1i3 = figure.text(-.0020,.015, '%2.0f' % float(freqi1i3array[discretemoment]), color='green', fontsize=4)
        else:
            self.bottomtexti1i3 = figure.text(-.0020,.015, '%f' % float(freqi1i3array[discretemoment]), color='red', fontsize=4)
        return self.bottomtexti1i3

    def makebottomtexthinge(figure, discretemoment, bottomtexthinge):
        if float(freqHingearray[discretemoment]) == 0:
            self.bottomtexthinge = figure.text(.0045,.015, '%2.0f' % float(freqHingearray[discretemoment]), color='black', fontsize=4)
        elif float(freqHingearray[discretemoment]) > 0:
            self.bottomtexthinge = figure.text(.0045,.015, '%2.0f' % float(freqHingearray[discretemoment]), color='green', fontsize=4)
        else:
            self.bottomtexthinge = figure.text(.0045,.015, '%f' % float(freqHingearray[discretemoment]), color='red', fontsize=4)
        return self.bottomtexthinge

    def maketimer(figure, discretemoment, timer):
        self.timer = figure.text(0, .02, 'Time: %3.2f seconds' % float(discretemoment/100.00), color='Black', fontsize=6, horizontalalignment = 'center')
        return self.timer

    def resetn3(n3graph):
        self.n3graph.remove()
        self.n3graph = figure.add_subplot(gs[13:14,1:8])
        #plt.tight_layout()
        return self.n3graph

    def resethinge(hingegraph):
        self.hingegraph.remove()
        self.hingegraph = figure.add_subplot(gs[17:18,1:8])
        #plt.tight_layout()
        return self.hingegraph

    def reseti2(i2graph):
        self.i2graph.remove()
        self.i2graph = figure.add_subplot(gs[11:12,1:8])
        #plt.tight_layout()
        return self.i2graph

    def reseti1i3(i1i3graph):
        self.i1i3graph.remove()
        self.i1i3graph = figure.add_subplot(gs[15:16,1:8])
        #plt.tight_layout()
        return self.i1i3graph

    ######################################################
    ##                 Create Video                     ##
    ######################################################
    def createMovie(self):
        #set which movie writer class is being used
        if animation.FFMpegWriter.isAvailable():
            self.movieWriterClass = animation.FFMpegWriter
        elif animation.AVConvWriter.isAvailable():
            self.movieWriterClass = animation.AVConvWriter
        else:
            sys.stderr.write('video converter missing')
            exit()

        #For each discretemoment/time step, save the figure as a frame in the .mp4 file
        self.moviewriter = self.movieWriterClass(fps = self.fps, bitrate = self.br)
        with self.moviewriter.saving(plt.gcf(), self.moviefile, self.dpi):
            i = 1 #first frame, must be 1<i<850
            frames = 850 #must be frames<=850
            while( i < frames):
                sys.stdout.flush()
                progress = str(round((i*100)/frames))
                sys.stdout.write(progress + '% complete \r' )
                #makefigures(i2graph, n3graph, i1i3graph, hingegraph, i)
                self.odontophore.remove()
                self.i1i3top.remove()
                self.i1i3bottom.remove()
                self.i1i3top1.remove()
                self.i1i3bottom1.remove()
                self.i1i3top2.remove()
                self.i1i3bottom2.remove()
                self.i1i3top3.remove()
                self.i1i3bottom3.remove()
                self.i1i3top4.remove()
                self.i1i3bottom4.remove()
                self.i2top.remove()
                self.i2bottom.remove()
                #self.vertlinei2.remove()
                #self.vertlinen3.remove()
                #self.vertlinei1i3.remove()
                #self.vertlinehinge.remove()
                #self.toptextI2.remove()
                #self.toptextN3.remove()
                #self.toptexti1i3.remove()
                #self.toptexthinge.remove()
                #self.bottomtexti2.remove()
                #self.bottomtextN3.remove()
                #self.bottomtexti1i3.remove()
                #self.bottomtexthinge.remove()
                #self.timer.remove()
                #self.toptextI2 = maketoptextI2(a, i, toptextI2)
                #self.toptextN3 = maketoptextN3(a, i, toptextN3)
                #self.toptexti1i3 = maketoptexti1i3(a, i, toptexti1i3)
                #self.toptexthinge = maketoptexthinge(a, i, toptexthinge)
                #self.bottomtexti2 = makebottomtexti2(a, i, bottomtexti2)
                #self.bottomtextN3 = makebottomtextN3(a, i, bottomtextN3)
                #self.bottomtexti1i3 = makebottomtexti1i3(a, i, bottomtexti1i3)
                #self.bottomtexthinge = makebottomtexthinge(a, i, bottomtexthinge)
                #self.timer = maketimer(a, i, timer)
                #self.vertlinei2 = makevertlinei2(i2graph, i)
                #self.vertlinen3 = makevertlinen3(n3graph, i)
                #self.vertlinei1i3 = makevertlinei1i3(i1i3graph, i)
                #self.vertlinehinge = makevertlinehinge(hingegraph, i)
                self.odontophore = self.createOdontophore(i, self.odontophore)
                self.i1i3top = self.createi1i3top(i, self.i1i3top)
                self.i1i3bottom = self.createi1i3bottom(i, self.i1i3bottom)
                self.i1i3top1 = self.createi1i3top1(i, self.i1i3top1)
                self.i1i3bottom1 = self.createi1i3bottom1(i, self.i1i3bottom1)
                self.i1i3top2 = self.createi1i3top2(i, self.i1i3top2)
                self.i1i3bottom2 = self.createi1i3bottom2(i, self.i1i3bottom2)
                self.i1i3top3 = self.createi1i3top3(i, self.i1i3top3)
                self.i1i3bottom3 = self.createi1i3bottom3(i, self.i1i3bottom3)
                self.i1i3top4 = self.createi1i3top4(i, self.i1i3top4)
                self.i1i3bottom4 = self.createi1i3bottom4(i, self.i1i3bottom4)
                self.i2top = self.createi2top(i, self.i2top)
                self.i2bottom = self.createi2bottom(i, self.i2bottom)
                self.a.add_artist(self.odontophore)
                self.a.add_artist(self.i1i3top)
                self.a.add_artist(self.i1i3bottom)
                self.a.add_artist(self.i1i3top1)
                self.a.add_artist(self.i1i3bottom1)
                self.a.add_artist(self.i1i3top2)
                self.a.add_artist(self.i1i3bottom2)
                self.a.add_artist(self.i1i3top3)
                self.a.add_artist(self.i1i3bottom3)
                self.a.add_artist(self.i1i3top4)
                self.a.add_artist(self.i1i3bottom4)
                self.a.add_artist(self.i2top)
                self.a.add_artist(self.i2bottom)
                self.moviewriter.grab_frame()
                i = i+5

#return mp4
myAnimator = Animator()
