import os
import csv
import matplotlib
import sys
sys.stdout.write('Starting... \r')
matplotlib.use('Agg') # force matplotlib to not use XWindows backend
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
#import plotly.plotly as py
#import plotly.tools as tls

#os.chdir("/Users/tatekeller/Documents/Aplysia Feeding Kinetic Model/AplysiaFeedingKineticModel")
#The following code reads the Izhikevich.txt file and converts it to a csv file
tabDelimitedFile = open("Izhikevich.txt", "r")
splitFileToList = "Izhikevich.txt".split(".")
listInput = []
for line in tabDelimitedFile:
    listLine = line.replace("\n","").split("\t")
    listInput.append(listLine)

tabDelimitedFile.close()
sys.stdout.write('File Read \r')
newFile = open(splitFileToList[0] + "-csvTranslated.csv", "w")
for line in listInput:
    writeLineForNewFile = ",".join(line)
    newFile.write(writeLineForNewFile + "\n")

newFile.close()
sys.stdout.write('New File Written \r')

with open("Izhikevich-csvTranslated.csv","r") as file:
    reader = csv.reader(file)
    data = list(reader)
    row_count = len(data)

#The following code reads the Izhikevich-csvTranslated.csv and saves the time in an array list "time"
time = []
i = 0
with open("Izhikevich-csvTranslated.csv", "r") as file:
     row = csv.reader(file)
     while i < row_count:
        for column in row:
            time.insert(i,column[0])
            i+=1
sys.stdout.write('Time Written \r')
#The following line prints the time array to confirm that it contains the correct values
#print(time)

#The following code reads the Izhikevich-csvTranslated.csv and saves v in an array list "membranePotential"
membranePotential = []
i = 0
with open("Izhikevich-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < row_count:
        for column in row:
            membranePotential.insert(i,column[1])
            i+=1
sys.stdout.write('Membrane Potential Written \r')
#The following line prints the time array to confirm that it contains the correct values
#print(membranePotential)

#The following code reads the Izhikevich-csvTranslated.csv and saves u in an array list "membraneRecovery"
membraneRecovery = []
i = 0
with open("Izhikevich-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < row_count:
        for column in row:
            membraneRecovery.insert(i,column[2])
            i+=1
sys.stdout.write('Membrane Recovery Written \r')
#The following code plots the entire figure
figure = plt.figure(figsize = (18,18))
timeticks = np.arange(0,8.5,1)

#The figure is a Grid with 2 rows, 1 columns.
gs = gridspec.GridSpec(2, 1)

#The following code subplots the membrane Potential vs time graph
potenialgraph = figure.add_subplot(gs[0,0])
potentialticks = np.arange(-90, 40, 10)
potenialgraph.plot(time[1:row_count],membranePotential[1:row_count])#first .008 seconds
potenialgraph.set_title('Membrane Potential')
potenialgraph.set_xlabel(time[0] + ' (milliseconds)')
potenialgraph.set_ylabel(membranePotential[0] + 'mV')
potenialgraph.axis([0, 8.5, -90, 40])
potenialgraph.set_xticks(timeticks)
potenialgraph.set_yticks(potentialticks)
potenialgraph.grid(True)

#The following code subplots the membrane recovery vs time graph
recoverygraph = figure.add_subplot(gs[1,0])
recoveryticks = np.arange(-20, 20, 2)
recoverygraph.plot(time[1:row_count],membraneRecovery[1:row_count])
recoverygraph.set_title('Membrane Recovery')
recoverygraph.set_xlabel(time[0] + ' (milliseconds)')
recoverygraph.set_ylabel(membraneRecovery[0])
recoverygraph.axis([0, 8.5, -20, 20])
recoverygraph.set_xticks(timeticks)
recoverygraph.set_yticks(recoveryticks)
recoverygraph.grid(True)

#Show the entire figure
#plt.show()

#Save Figure to PDF
figure.tight_layout()
figure.savefig("plot.pdf")
