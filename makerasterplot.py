import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg') # force matplotlib to not use XWindows backend
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# read in the data
with open("rasterplotinfo.csv","r") as file:
    data = pd.read_csv(file)

# find spike times and compute average firing frequencies
spike_times = {}
firing_freq = {}
for neuron in data.columns[1:]:
    spike_times[neuron] = data['time'][data[neuron] == 1].values
    firing_freq[neuron] = (spike_times[neuron].size - 1) / np.ptp(spike_times[neuron])

gs = gridspec.GridSpec(2, 1)
figure = plt.figure(figsize = (18,18))
rasterplot = figure.add_subplot(gs[1,0])

rasterplot.eventplot(
    (
        (data['time'][data['B7'] == 1]),
        (data['time'][data['B43'] == 1]),
        (data['time'][data['B10'] == 1]),
        (data['time'][data['B38'] == 1]),
        (data['time'][data['B9'] == 1]),
        (data['time'][data['B6'] == 1]),
        (data['time'][data['B3'] == 1]),
        (data['time'][data['B8b'] == 1]),
        (data['time'][data['B8a'] == 1]),
        (data['time'][data['B61/62'] == 1]),
        (data['time'][data['B31/32'] == 1]),
    ),
    linelengths = [.9,.9,.9,.9,.9,.9,.9,.9,.9,.9,.9],
    linewidths = [2,2,2,2,2,2,2,2,2,2,2],
    colors = ['hotpink', 'gray', 'orangered', 'purple', 'dodgerblue', 'dodgerblue', 'gold', 'mediumseagreen', 'mediumseagreen', 'dimgray', 'dimgray'],
    )

timeticks = np.arange(0,8,.5)
neuralticks = np.arange(0,11,1)
rasterplot.set_yticks(neuralticks)
rasterplot.set_yticklabels([])
rasterplot.set_xticks(timeticks)
rasterplot.tick_params(labelsize = 20)
rasterplot.text(3.75,-3, 'Time',color = 'black', fontsize=20,horizontalalignment = 'center',fontweight='bold')
rasterplot.grid(True)
rasterplot.text(-0.2,10.85, 'Neuron',color = 'black', fontsize=20,horizontalalignment = 'right',fontweight='bold')
rasterplot.text(-0.2,9.85, 'B31/32',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.2,8.85, 'B61/62',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.2,7.85, 'B8a',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.2,6.85, 'B8b',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.2,5.85, 'B3',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.2,4.85, 'B6',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.2,3.85, 'B9',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.2,2.85, 'B38',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.2,1.85, 'B10',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.2,0.85, 'B43',color = 'black', fontsize=20,horizontalalignment = 'right')
rasterplot.text(-0.2,-0.15, 'B7',color = 'black', fontsize=20,horizontalalignment = 'right')

rasterplot.text(7.55,10.85, 'Frequency',color = 'black', fontsize=20,horizontalalignment = 'left',fontweight='bold')
rasterplot.text(7.55,9.85, '%2.1f Hz' % (firing_freq['B31/32']),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(7.55,8.85, '%2.1f Hz' % (firing_freq['B61/62']),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(7.55,7.85, '%2.1f Hz' % (firing_freq['B8a']),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(7.55,6.85, '%2.1f Hz' % (firing_freq['B8b']),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(7.55,5.85, '%2.1f Hz' % (firing_freq['B3']),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(7.55,4.85, '%2.1f Hz' % (firing_freq['B6']),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(7.55,3.85, '%2.1f Hz' % (firing_freq['B9']),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(7.55,2.85, '%2.1f Hz' % (firing_freq['B38']),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(7.55,1.85, '%2.1f Hz' % (firing_freq['B10']),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(7.55,0.85, '%2.1f Hz' % (firing_freq['B43']),color = 'black', fontsize=20,horizontalalignment = 'left')
rasterplot.text(7.55,-0.15, '%2.1f Hz' % (firing_freq['B7']),color = 'black', fontsize=20,horizontalalignment = 'left')

plt.title('Raster Plot', fontsize=40)
figure.savefig("RasterPlot.pdf")
