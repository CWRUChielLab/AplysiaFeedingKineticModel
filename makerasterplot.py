import os
import csv
import sys
sys.stdout.write('Starting... \r')
import matplotlib
#matplotlib.use('Agg') # force matplotlib to not use XWindows backend
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

#The following code reads the rasterplotinfo.txt file and converts it to a csv file
tabDelimitedFile = open("rasterplotinfo.txt", "r")
splitFileToList = "rasterplotinfo.txt".split(".")
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

with open("rasterplotinfo-csvTranslated.csv","r") as file:
    reader = csv.reader(file)
    data = list(reader)
    row_count = len(data)

#The following code reads the rasterplotinfo-csvTranslated.csv and saves the time in an array list "time"
time = []
i = 0
with open("rasterplotinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < row_count:
        for column in row:
            time.insert(i,column[0])
            i+=1
sys.stdout.write('Time Saved \r')
# [0]- B31/32, [1] - B61/62, [2] - B8a,b, [3] - B3, [4] - B6, [5] - B9, [6] - B38, [7] - B10, [8] - B43, [9] - B7

#The following code reads the rasterplotinfo-csvTranslated.csv and saves the time in an array list "time"
b31 = []
i = 0
with open("rasterplotinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < row_count:
        for column in row:
            b31.insert(i,column[1])
            i+=1
sys.stdout.write('B31/32 Saved \r')
#The following code reads the rasterplotinfo-csvTranslated.csv and saves the time in an array list "time"
b61 = []
i = 0
with open("rasterplotinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < row_count:
        for column in row:
            b61.insert(i,column[2])
            i+=1
sys.stdout.write('B61/62 Saved \r')
#The following code reads the rasterplotinfo-csvTranslated.csv and saves the time in an array list "time"
b8 = []
i = 0
with open("rasterplotinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < row_count:
        for column in row:
            b8.insert(i,column[3])
            i+=1
sys.stdout.write('B8 Saved \r')
#The following code reads the rasterplotinfo-csvTranslated.csv and saves the time in an array list "time"
b3 = []
i = 0
with open("rasterplotinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < row_count:
        for column in row:
            b3.insert(i,column[4])
            i+=1
sys.stdout.write('B3 Saved \r')
#The following code reads the rasterplotinfo-csvTranslated.csv and saves the time in an array list "time"
b6 = []
i = 0
with open("rasterplotinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < row_count:
        for column in row:
            b6.insert(i,column[5])
            i+=1
sys.stdout.write('B6 Saved \r')
#The following code reads the rasterplotinfo-csvTranslated.csv and saves the time in an array list "time"
b9 = []
i = 0
with open("rasterplotinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < row_count:
        for column in row:
            b9.insert(i,column[6])
            i+=1
sys.stdout.write('B9 Saved \r')
#The following code reads the rasterplotinfo-csvTranslated.csv and saves the time in an array list "time"
b38 = []
i = 0
with open("rasterplotinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < row_count:
        for column in row:
            b38.insert(i,column[7])
            i+=1
sys.stdout.write('B38 Saved \r')
#The following code reads the rasterplotinfo-csvTranslated.csv and saves the time in an array list "time"
b10 = []
i = 0
with open("rasterplotinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < row_count:
        for column in row:
            b10.insert(i,column[8])
            i+=1
sys.stdout.write('B10 Saved \r')
#The following code reads the rasterplotinfo-csvTranslated.csv and saves the time in an array list "time"
b43 = []
i = 0
with open("rasterplotinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < row_count:
        for column in row:
            b43.insert(i,column[9])
            i+=1
sys.stdout.write('B43 Saved \r')
#The following code reads the rasterplotinfo-csvTranslated.csv and saves the time in an array list "time"
b7 = []
i = 0
with open("rasterplotinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < row_count:
        for column in row:
            b7.insert(i,column[10])
            i+=1
sys.stdout.write('B7 Saved \r')

#The following code reads the rasterplotinfo-csvTranslated.csv and saves the time in an array list "time"
b8b = []
i = 0
with open("rasterplotinfo-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < row_count:
        for column in row:
            b8b.insert(i,column[11])
            i+=1
sys.stdout.write('B8b Saved \r')

b31moment = []
for i in range(1,row_count):
    if(float(b31[i]) == 1):
        b31moment.append(time[i])

b61moment = []
for i in range(1,row_count):
    if(float(b61[i]) == 1):
        b61moment.append(time[i])

b8moment = []
for i in range(1,row_count):
    if(float(b8[i]) == 1):
        b8moment.append(time[i])

b3moment = []
for i in range(1,row_count):
    if(float(b3[i]) == 1):
        b3moment.append(time[i])

b6moment = []
for i in range(1,row_count):
    if(float(b6[i]) == 1):
        b6moment.append(time[i])

b9moment = []
for i in range(1,row_count):
    if(float(b9[i]) == 1):
        b9moment.append(time[i])

b38moment = []
for i in range(1,row_count):
    if(float(b38[i]) == 1):
        b38moment.append(time[i])

b10moment = []
for i in range(1,row_count):
    if(float(b10[i]) == 1):
        b10moment.append(time[i])

b43moment = []
for i in range(1,row_count):
    if(float(b43[i]) == 1):
        b43moment.append(time[i])

b7moment = []
for i in range(1,row_count):
    if(float(b7[i]) == 1):
        b7moment.append(time[i])

b8bmoment = []
for i in range(1,row_count):
    if(float(b8b[i]) == 1):
        b8bmoment.append(time[i])

gs = gridspec.GridSpec(1, 1)
figure = plt.figure(figsize = (18,18))
rasterplot = figure.add_subplot(gs[0,0])

x = np.array(b31moment[:])
b31momentarray = x.astype(np.float)

x = np.array(b61moment[:])
b61momentarray = x.astype(np.float)

x = np.array(b8moment[:])
b8momentarray = x.astype(np.float)

x = np.array(b3moment[:])
b3momentarray = x.astype(np.float)

x = np.array(b6moment[:])
b6momentarray = x.astype(np.float)

x = np.array(b9moment[:])
b9momentarray = x.astype(np.float)

x = np.array(b38moment[:])
b38momentarray = x.astype(np.float)

x = np.array(b10moment[:])
b10momentarray = x.astype(np.float)

x = np.array(b43moment[:])
b43momentarray = x.astype(np.float)

x = np.array(b7moment[:])
b7momentarray = x.astype(np.float)

x = np.array(b8bmoment[:])
b8bmomentarray = x.astype(np.float)

rasterplot.eventplot(((b7momentarray[:]), (b43momentarray[:]), (b10momentarray[:]),(b38momentarray[:]), (b9momentarray[:]), (b6momentarray[:]),(b3momentarray[:]),(b8bmomentarray[:]),(b8momentarray[:]),(b61momentarray[:]),(b31momentarray[:])), linelengths = [.9,.9,.9,.9,.9,.9,.9,.9,.9,.9,.9], colors = ['dimgray', 'dodgerblue', 'dodgerblue', 'dodgerblue', 'dodgerblue', 'dodgerblue', 'dodgerblue', 'coral', 'coral', 'mediumseagreen', 'mediumseagreen'], linewidths = [.7,.7,.7,.7,.7,.7,.7,.7,.7,.7,.7])

timeticks = np.arange(0,9,.5)
neuralticks = np.arange(0,11,1)
rasterplot.set_yticks(neuralticks)
rasterplot.set_yticklabels([])
rasterplot.set_xticks(timeticks)
rasterplot.grid(True)
rasterplot.text(-0.5,9.95, 'B31/32',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.5,8.95, 'B61/62',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.5,7.95, 'B8a',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.5,6.95, 'B8b',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.5,5.95, 'B3',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.5,4.95, 'B6',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.5,3.95, 'B9',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.5,2.95, 'B38',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.5,1.95, 'B10',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.5,0.95, 'B43',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.5,-0.05, 'B7',color = 'black', fontsize=20,horizontalalignment = 'right')

rasterplot.text(8.55,9.95, '%2.1f Hz' % ((len(b31momentarray)-1)/(b31momentarray[-1]-b31momentarray[0])),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(8.55,8.95, '%2.1f Hz' % ((len(b61momentarray)-1)/(b61momentarray[-1]-b61momentarray[0])),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(8.55,7.95, '%2.1f Hz' % ((len(b8momentarray)-1)/(b8momentarray[-1]-b8momentarray[0])),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(8.55,6.95, '%2.1f Hz' % ((len(b8bmomentarray)-1)/(b8bmomentarray[-1]-b8bmomentarray[0])),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(8.55,5.95, '%2.1f Hz' % ((len(b3momentarray)-1)/(b3momentarray[-1]-b3momentarray[0])),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(8.55,4.95, '%2.1f Hz' % ((len(b6momentarray)-1)/(b6momentarray[-1]-b6momentarray[0])),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(8.55,3.95, '%2.1f Hz' % ((len(b9momentarray)-1)/(b9momentarray[-1]-b9momentarray[0])),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(8.55,2.95, '%2.1f Hz' % ((len(b38momentarray)-1)/(b38momentarray[-1]-b38momentarray[0])),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(8.55,1.95, '%2.1f Hz' % ((len(b10momentarray)-1)/(b10momentarray[-1]-b10momentarray[0])),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(8.55,0.95, '%2.1f Hz' % ((len(b43momentarray)-1)/(b43momentarray[-1]-b43momentarray[0])),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(8.55,-0.05, '%2.1f Hz' % ((len(b7momentarray)-1)/(b7momentarray[-1]-b7momentarray[0])),color = 'black', fontsize=20,horizontalalignment = 'left')





plt.title('Raster Plot', fontsize=40)
figure.savefig("RasterPlot.pdf")