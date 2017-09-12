import os
import csv
import matplotlib.pyplot as plt
#os.chdir("/Users/tatekeller/Documents/Aplysia Feeding Kinetic Model/AplysiaFeedingKineticModel")
#The following code reads the SlugOutput2.txt file and converts it to a csv file
tabDelimitedFile = open("SlugOutput2.txt", "r")
splitFileToList = "SlugOutput2.txt".split(".")
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
#The following code reads the SlugOutput2-csvTranslated.csv and saves the time in an array list "time"
time = []
i = 0
with open("SlugOutput2-csvTranslated.csv", "r") as file:
     row = csv.reader(file)
     while i < 850:
        for column in row:
            time.insert(i,column[0])
            i+=1

#The following line prints the time array to confirm that it contains the correct values
print(time)

#The following code reads the SlugOutput2-csvTranslated.csv and saves the position in an array list "position"
position = []
i = 0
with open("SlugOutput2-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            position.insert(i,column[1])
            i+=1

#The following line prints the time array to confirm that it contains the correct values
print(position)

#The following code reads the SlugOutput2-csvTranslated.csv and saves the radius in an array list "radius"
radius = []
i = 0
with open("SlugOutput2-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            radius.insert(i,column[2])
            i+=1

#The following line prints the time array to confirm that it contains the correct values
print(radius)

#The following code reads the SlugOutput2-csvTranslated.csv and saves the angle in an array list "angle"
angle = []
i = 0
with open("SlugOutput2-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            angle.insert(i,column[3])
            i+=1

#The following line prints the time array to confirm that it contains the correct values
print(angle)

#The following code reads the SlugOutput2-csvTranslated.csv and saves the seaweedforce in an array list "seaweedforce"
seaweedforce = []
i = 0
with open("SlugOutput2-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            seaweedforce.insert(i,column[14])
            i+=1

#The following line prints the time array to confirm that it contains the correct values
print(seaweedforce)

#The following code reads the SlugOutput2-csvTranslated.csv and saves the Odontophore input in an array list "freqN3"
freqN3 = []
i = 0
with open("SlugOutput2-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            freqN3.insert(i,column[8])
            i+=1

#The following line prints the time array to confirm that it contains the correct values
print(freqN3)

#The following code reads the SlugOutput2-csvTranslated.csv and saves the Hinge input in an array list "freqHinge"
freqHinge = []
i = 0
with open("SlugOutput2-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            freqHinge.insert(i,column[9])
            i+=1

#The following line prints the time array to confirm that it contains the correct values
print(freqHinge)

#The following code reads the SlugOutput2-csvTranslated.csv and saves the I2 input in an array list "freqI2"
freqI2 = []
i = 0
with open("SlugOutput2-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            freqI2.insert(i,column[6])
            i+=1

#The following line prints the time array to confirm that it contains the correct values
print(freqI2)

#The following code reads the SlugOutput2-csvTranslated.csv and saves the I1/I3 input in an array list "freqI1I3"
freqI1I3 = []
i = 0
with open("SlugOutput2-csvTranslated.csv", "r") as file:
    row = csv.reader(file)
    while i < 850:
        for column in row:
            freqI1I3.insert(i,column[7])
            i+=1

#The following line prints the time array to confirm that it contains the correct values
print(freqI1I3)

#The following code plots the entire figure
plt.figure(1)

#The following code subplots the position vs time graph
plt.subplot(811)
plt.plot(time[1:850],position[1:850])
plt.xlabel(time[0])
plt.ylabel(position[0])
plt.axis([0, 9, -0.006, 0.004])
plt.grid(True)
#The following code subplots the Odontophore major axis vs time graph
plt.subplot(812)
plt.plot(time[1:850],radius[1:850])
plt.xlabel(time[0])
plt.ylabel(radius[0])
plt.axis([0, 9, 0, 0.009])
plt.grid(True)
#The following code subplots the Angle vs time graph
plt.subplot(813)
plt.plot(time[1:850],angle[1:850])
plt.xlabel(time[0])
plt.ylabel(angle[0])
plt.axis([0, 9, -10, 100])
plt.grid(True)
#The following code subplots the Seaweed Force vs time graph
plt.subplot(814)
plt.plot(time[1:850],seaweedforce[1:850])
plt.xlabel(time[0])
plt.ylabel(seaweedforce[0])
plt.axis([0, 9, -20, 120])
plt.grid(True)
#The following code subplots the Odontophore input vs time graph
plt.subplot(815)
plt.plot(time[1:850],freqN3[1:850])
plt.xlabel(time[0])
plt.ylabel(freqN3[0])
plt.axis([0, 9, -5, 35])
plt.grid(True)
#The following code subplots the Hinge input vs time graph
plt.subplot(816)
plt.plot(time[1:850],freqHinge[1:850])
plt.xlabel(time[0])
plt.ylabel(freqHinge[0])
plt.axis([0, 9, -5, 25])
plt.grid(True)
#The following code subplots the I2 input vs time graph
plt.subplot(817)
plt.plot(time[1:850],freqI2[1:850])
plt.xlabel(time[0])
plt.ylabel(freqI2[0])
plt.axis([0, 9, -5, 25])
plt.grid(True)
#The following code subplots the I1/I3 vs time graph
plt.subplot(818)
plt.plot(time[1:850],freqI1I3[1:850])
plt.xlabel(time[0])
plt.ylabel(freqI1I3[0])
plt.axis([0, 9, -5, 25])
plt.grid(True)
#Show the entire figure
plt.show()
