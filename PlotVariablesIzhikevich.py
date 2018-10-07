import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg') # force matplotlib to not use XWindows backend
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# read in the data
with open("Izhikevich.csv","r") as file:
    data = pd.read_csv(file)

#The following code plots the entire figure
figure = plt.figure(figsize = (18,18))
timeticks = np.arange(0,8.5,1)

#The figure is a Grid with 2 rows, 1 columns.
gs = gridspec.GridSpec(2, 1)

#The following code subplots the membrane Potential vs time graph
potenialgraph = figure.add_subplot(gs[0,0])
potentialticks = np.arange(-90, 40, 10)
potenialgraph.plot(data['time'], data['MembranePotentialo'])
potenialgraph.set_title('Membrane Potential')
potenialgraph.set_xlabel('Time (milliseconds)')
potenialgraph.set_ylabel('Membrane Potential (mV)')
potenialgraph.axis([0, 8.5, -90, 40])
potenialgraph.set_xticks(timeticks)
potenialgraph.set_yticks(potentialticks)
potenialgraph.grid(True)

#The following code subplots the membrane recovery vs time graph
recoverygraph = figure.add_subplot(gs[1,0])
recoveryticks = np.arange(-20, 20, 2)
recoverygraph.plot(data['time'], data['MembraneRecoveryo'])
recoverygraph.set_title('Membrane Recovery')
recoverygraph.set_xlabel('Time (milliseconds)')
recoverygraph.set_ylabel('Membrane Recovery')
recoverygraph.axis([0, 8.5, -20, 20])
recoverygraph.set_xticks(timeticks)
recoverygraph.set_yticks(recoveryticks)
recoverygraph.grid(True)

#Show the entire figure
#plt.show()

#Save Figure to PDF
figure.tight_layout()
figure.savefig("plot.pdf")
