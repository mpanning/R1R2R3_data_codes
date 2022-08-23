"""
Back azimuth figure
"""

from obspy.core import read, UTCDateTime
from obspy.signal.filter import envelope
from obspy.signal.rotate import rotate_ne_rt
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert
from matplotlib.patches import Rectangle

datafile = './data/S1222a.BH.vbbrotate.deglitched.mseed'
origin = UTCDateTime(2022,5,4,23,23,6,570000)
reftime = UTCDateTime(2022,5,4,23,27,45,690000)

cfreqs = [0.021, 0.025, 0.03, 0.036, 0.043,
          0.05, 0.06, 0.072, 0.086]
cfreqs = np.array(cfreqs)
baz_angles = np.arange(0, 360, 2)
hwidth = 0.3

# Set up plot grid
fig = plt.figure(figsize=(8, 9))
grid = plt.GridSpec(9, 3, hspace=0, wspace=0.2)
axs = [plt.subplot(grid[:3, 0]), plt.subplot(grid[:3, 1]),
       plt.subplot(grid[4:7, 0]), plt.subplot(grid[4:7, 1]),
       plt.subplot(grid[8, :2]), plt.subplot(grid[:7,2])]
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
labelxoff = 0.03
labelyoff = 0.9

# R1 window
ax = axs[0]
plt.sca(ax)

SWpicks = [490.25, 490.84, 501.55, 512.42, 521.27, 529.21,
           536.00, 585.77, 586.76]

corr_angles = np.zeros((cfreqs.size, baz_angles.size))
st = read(datafile)

for i, cfreq in enumerate(cfreqs):
    print('Working on frequency {}'.format(cfreq))
    stf = st.copy()
    stf.filter('bandpass', freqmin=cfreq*(1-hwidth),
               freqmax=cfreq*(1+hwidth), corners=2, zerophase=True)
    stime = reftime + SWpicks[i] - 2./cfreq
    etime = stime + 4./cfreq
    stf.trim(starttime=stime, endtime=etime)
    dataz_orig = stf.select(channel='BHZ')[0].data
    dataz_hilbert = np.imag(hilbert(dataz_orig))
    # baz_angles = np.arange(0, 360, 5)
    datan = stf.select(channel='BHN')[0].data
    datae = stf.select(channel='BHE')[0].data
    for j, baz in enumerate(baz_angles):
        print('Working on baz {}'.format(baz))
        datar, datat = rotate_ne_rt(datan, datae, baz)
        corr_angles[i,j] = -1.*np.corrcoef(dataz_hilbert, datar)[0, 1]

    plt.plot(baz_angles, corr_angles[i,:], color=colors[i],
             label='{} Hz'.format(cfreq), linestyle='dashed')

corr_ave = np.mean(corr_angles, axis=0)
corr_ave_high = np.mean(corr_angles[1:8,:], axis=0)

stf = st.copy()
stf.filter('bandpass', freqmin=0.018, freqmax=0.1, corners=2, zerophase=True)
stime = reftime + 450.0
etime = reftime + 650.0
stf.trim(starttime=stime, endtime=etime)
dataz_orig = stf.select(channel='BHZ')[0].data
dataz_hilbert = np.imag(hilbert(dataz_orig))
datan = stf.select(channel='BHN')[0].data
datae = stf.select(channel='BHE')[0].data
corr_bb = np.zeros(baz_angles.size)
for i, baz in enumerate(baz_angles):
    datar, datat = rotate_ne_rt(datan, datae, baz)
    corr_bb[i] = -1.*np.corrcoef(dataz_hilbert, datar)[0, 1]
plt.plot(baz_angles, corr_ave, 'k-', label='Mean')
plt.plot(baz_angles, corr_ave_high, 'r-', label='High Qual')
plt.plot(baz_angles, corr_bb, 'b-', label='Broadband')
# plt.legend()

#Add in MQS baz pick
MQSbaz = 109
MQSbaz_sigma = 7.44
plt.plot([MQSbaz, MQSbaz], [-1, 1], 'k-.', label="MQS baz")
ax.add_patch(Rectangle((MQSbaz - 2.*MQSbaz_sigma, -1), 4.*MQSbaz_sigma, 2,
                       color='grey', alpha=0.5))

plt.ylim(-1,1)
plt.xlim(0,360)
ax.set_title("R1 window")
ax.set_ylabel("Correlation")
ax.text(labelxoff, labelyoff, "A", transform=ax.transAxes, fontsize=16)

# Plot legend
h, l = ax.get_legend_handles_labels()
axs[4].legend(h,l, borderaxespad=0, loc=10, ncol=3)
axs[4].axis("off")

# R2 window
ax = axs[1]
plt.sca(ax)

SWpicks = [6375.26, 6373.66, 6375.01, 6376.71, 6377.97, 6383.81]
st = read(datafile)

for i, cfreq in enumerate(cfreqs[:6]):
    print('Working on frequency {}'.format(cfreq))
    stf = st.copy()
    stf.filter('bandpass', freqmin=cfreq*(1-hwidth),
               freqmax=cfreq*(1+hwidth), corners=2, zerophase=True)
    stime = reftime + SWpicks[i] - 2./cfreq
    etime = stime + 4./cfreq
    stf.trim(starttime=stime, endtime=etime)
    dataz_orig = stf.select(channel='BHZ')[0].data
    dataz_hilbert = np.imag(hilbert(dataz_orig))
    # baz_angles = np.arange(0, 360, 5)
    datan = stf.select(channel='BHN')[0].data
    datae = stf.select(channel='BHE')[0].data
    for j, baz in enumerate(baz_angles):
        print('Working on baz {}'.format(baz))
        datar, datat = rotate_ne_rt(datan, datae, baz)
        corr_angles[i,j] = -1.*np.corrcoef(dataz_hilbert, datar)[0, 1]

    plt.plot(baz_angles, corr_angles[i,:], color=colors[i],
             label='{} Hz'.format(cfreq), linestyle='dashed')

corr_ave = np.mean(corr_angles[:6,:], axis=0)
corr_ave_high = np.mean(corr_angles[1:5,:], axis=0)

stf = st.copy()
stf.filter('bandpass', freqmin=0.021, freqmax=0.05, corners=2, zerophase=True)
stime = reftime + 6250.0
etime = reftime + 6450.0
stf.trim(starttime=stime, endtime=etime)
dataz_orig = stf.select(channel='BHZ')[0].data
dataz_hilbert = np.imag(hilbert(dataz_orig))
datan = stf.select(channel='BHN')[0].data
datae = stf.select(channel='BHE')[0].data
corr_bb = np.zeros(baz_angles.size)
for i, baz in enumerate(baz_angles):
    datar, datat = rotate_ne_rt(datan, datae, baz)
    corr_bb[i] = -1.*np.corrcoef(dataz_hilbert, datar)[0, 1]
plt.plot(baz_angles, corr_ave, 'k-', label='Mean')
plt.plot(baz_angles, corr_ave_high, 'r-', label='High Qual')
plt.plot(baz_angles, corr_bb, 'b-', label='Broadband')
#Add in MQS baz pick
plt.plot([MQSbaz+180, MQSbaz+180], [-1, 1], 'k-.')
ax.add_patch(Rectangle((MQSbaz+180 - 2.*MQSbaz_sigma, -1), 4.*MQSbaz_sigma, 2,
                       color='grey', alpha=0.5))
plt.ylim(-1,1)
plt.xlim(0,360)
ax.text(labelxoff, labelyoff, "B", transform=ax.transAxes, fontsize=16)
ax.set_title("R2 window")
ax.yaxis.set_ticklabels([])

# R3 window
ax = axs[2]
plt.sca(ax)

SWpicks = [7872.31, 7863.56, 7864.21, 7882.36, 7882.72]
st = read(datafile)

for i, cfreq in enumerate(cfreqs[:5]):
    print('Working on frequency {}'.format(cfreq))
    stf = st.copy()
    stf.filter('bandpass', freqmin=cfreq*(1-hwidth),
               freqmax=cfreq*(1+hwidth), corners=2, zerophase=True)
    stime = reftime + SWpicks[i] - 2./cfreq
    etime = stime + 4./cfreq
    stf.trim(starttime=stime, endtime=etime)
    dataz_orig = stf.select(channel='BHZ')[0].data
    dataz_hilbert = np.imag(hilbert(dataz_orig))
    # baz_angles = np.arange(0, 360, 5)
    datan = stf.select(channel='BHN')[0].data
    datae = stf.select(channel='BHE')[0].data
    for j, baz in enumerate(baz_angles):
        print('Working on baz {}'.format(baz))
        datar, datat = rotate_ne_rt(datan, datae, baz)
        corr_angles[i,j] = -1.*np.corrcoef(dataz_hilbert, datar)[0, 1]

    plt.plot(baz_angles, corr_angles[i,:], color=colors[i],
             label='{} Hz'.format(cfreq), linestyle='dashed')

corr_ave = np.mean(corr_angles[:5,:], axis=0)
corr_ave_high = np.mean(corr_angles[3:5,:], axis=0)

stf = st.copy()
stf.filter('bandpass', freqmin=0.021, freqmax=0.05, corners=2, zerophase=True)
stime = reftime + 7825.0
etime = reftime + 7975.0
stf.trim(starttime=stime, endtime=etime)
dataz_orig = stf.select(channel='BHZ')[0].data
dataz_hilbert = np.imag(hilbert(dataz_orig))
datan = stf.select(channel='BHN')[0].data
datae = stf.select(channel='BHE')[0].data
corr_bb = np.zeros(baz_angles.size)
for i, baz in enumerate(baz_angles):
    datar, datat = rotate_ne_rt(datan, datae, baz)
    corr_bb[i] = -1.*np.corrcoef(dataz_hilbert, datar)[0, 1]
plt.plot(baz_angles, corr_ave, 'k-', label='Mean')
plt.plot(baz_angles, corr_ave_high, 'r-', label='High Qual')
plt.plot(baz_angles, corr_bb, 'b-', label='Broadband')
#Add in MQS baz pick
plt.plot([MQSbaz, MQSbaz], [-1, 1], 'k-.')
ax.add_patch(Rectangle((MQSbaz - 2.*MQSbaz_sigma, -1), 4.*MQSbaz_sigma, 2,
                       color='grey', alpha=0.5))
plt.ylim(-1,1)
plt.xlim(0,360)
ax.text(labelxoff, labelyoff, "C", transform=ax.transAxes, fontsize=16)
ax.set_title("R3 window")
ax.set_ylabel("Correlation")
ax.set_xlabel("Backazimuth")

# R2 window
ax = axs[3]
plt.sca(ax)

SWpicks = [6661.96, 6653.51, 6688.06, 6697.87, 6708.01, 6726.36]
st = read(datafile)

for i, cfreq in enumerate(cfreqs[:6]):
    print('Working on frequency {}'.format(cfreq))
    stf = st.copy()
    stf.filter('bandpass', freqmin=cfreq*(1-hwidth),
               freqmax=cfreq*(1+hwidth), corners=2, zerophase=True)
    stime = reftime + SWpicks[i] - 2./cfreq
    etime = stime + 4./cfreq
    stf.trim(starttime=stime, endtime=etime)
    dataz_orig = stf.select(channel='BHZ')[0].data
    dataz_hilbert = np.imag(hilbert(dataz_orig))
    # baz_angles = np.arange(0, 360, 5)
    datan = stf.select(channel='BHN')[0].data
    datae = stf.select(channel='BHE')[0].data
    for j, baz in enumerate(baz_angles):
        print('Working on baz {}'.format(baz))
        datar, datat = rotate_ne_rt(datan, datae, baz)
        corr_angles[i,j] = -1.*np.corrcoef(dataz_hilbert, datar)[0, 1]

    plt.plot(baz_angles, corr_angles[i,:], color=colors[i],
             label='{} Hz'.format(cfreq), linestyle='dashed')

corr_ave = np.mean(corr_angles[:6,:], axis=0)
corr_ave_high = np.mean(corr_angles[1:5,:], axis=0)

stf = st.copy()
stf.filter('bandpass', freqmin=0.021, freqmax=0.05, corners=2, zerophase=True)
stime = reftime + 6550.0
etime = reftime + 6750.0
stf.trim(starttime=stime, endtime=etime)
dataz_orig = stf.select(channel='BHZ')[0].data
dataz_hilbert = np.imag(hilbert(dataz_orig))
datan = stf.select(channel='BHN')[0].data
datae = stf.select(channel='BHE')[0].data
corr_bb = np.zeros(baz_angles.size)
for i, baz in enumerate(baz_angles):
    datar, datat = rotate_ne_rt(datan, datae, baz)
    corr_bb[i] = -1.*np.corrcoef(dataz_hilbert, datar)[0, 1]
plt.plot(baz_angles, corr_ave, 'k-', label='Mean')
plt.plot(baz_angles, corr_ave_high, 'r-', label='High Qual')
plt.plot(baz_angles, corr_bb, 'b-', label='Broadband')
#Add in MQS baz pick
plt.plot([MQSbaz+180, MQSbaz+180], [-1, 1], 'k-.')
ax.add_patch(Rectangle((MQSbaz+180 - 2.*MQSbaz_sigma, -1), 4.*MQSbaz_sigma, 2,
                       color='grey', alpha=0.5))
plt.ylim(-1,1)
plt.xlim(0,360)
ax.text(labelxoff, labelyoff, "D", transform=ax.transAxes, fontsize=16)
ax.set_title("Alternate R2")
ax.set_xlabel("Backazimuth")
ax.yaxis.set_ticklabels([])

# Plot up different baz estimates from different methods and different phases
ax = axs[5]
plt.sca(ax)

# R1 picks
# 3 JPL picks
plt.plot([108], [-0], label='R1 Ia', color=colors[0], marker="o",
         markersize=5)
plt.plot([68, 176], [-0, -0], color=colors[0])
plt.plot([146], [-0.1], label='R1 Ib', color=colors[0], marker="o",
         markersize=5)
plt.plot([94, 200], [-0.1, -0.1], color=colors[0], linestyle="dashed")
plt.plot([140], [-0.2], label='R1 Ic', color=colors[0], marker="o",
         markersize=5)
plt.plot([72, 202], [-0.2, -0.2], color=colors[0], linestyle="dotted")
# Carrasco
plt.plot([130], [-0.3], label='R1 II', color=colors[1], marker="o",
         markersize=5)
plt.plot([114, 137], [-0.3, -0.3], color=colors[1])
# DK
plt.plot([124.5], [-0.4], label='R1 III', color=colors[2], marker='o',
         markersize=5)
plt.plot([114.5, 139.7], [-0.4, -0.4], color=colors[2])
ax.text(270, -0.2, "R1", fontsize=12)
#Add in MQS baz pick
plt.plot([MQSbaz, MQSbaz], [0.1, -0.5], 'k-.')
ax.add_patch(Rectangle((MQSbaz - 2.*MQSbaz_sigma, -0.5), 4.*MQSbaz_sigma,
                       0.6, color='grey', alpha=0.5))

# R2 picks
plt.plot([0, 360], [-0.5, -0.5], color='k')
# 3 JPL picks
plt.plot([348], [-0.6], label='R2 Ia', color=colors[0], marker="s",
         markersize=5)
plt.plot([318, 360], [-0.6, -0.6], color=colors[0])
plt.plot([0, 18], [-0.6, -0.6], color=colors[0])
plt.plot([340], [-0.7], label='R2 Ib', color=colors[0], marker="s",
         markersize=5)
plt.plot([318, 360], [-0.7, -0.7], color=colors[0], linestyle="dashed")
plt.plot([0, 12], [-0.7, -0.7], color=colors[0], linestyle="dashed")
plt.plot([340], [-0.8], label='R2 Ic', color=colors[0], marker="s",
         markersize=5)
plt.plot([308, 360], [-0.8, -0.8], color=colors[0], linestyle="dotted")
plt.plot([0, 12], [-0.8, -0.8], color=colors[0], linestyle="dotted")
# Carrasco
plt.plot([293], [-0.9], label='R2 II', color=colors[1], marker="s",
         markersize=5)
plt.plot([266, 301], [-0.9, -0.9], color=colors[1])
# DK
plt.plot([197.1], [-1.0], label='R2 III', color=colors[2], marker='o',
         markersize=5)
plt.plot([188, 224], [-1.0, -1.0], color=colors[2])
ax.text(90, -0.8, "R2", fontsize=12)
#Add in MQS baz pick
plt.plot([MQSbaz+180, MQSbaz+180], [-0.5, -1.1], 'k-.')
ax.add_patch(Rectangle((MQSbaz+180 - 2.*MQSbaz_sigma, -1.1), 4.*MQSbaz_sigma,
                       0.6, color='grey', alpha=0.5))

# R3 picks
plt.plot([0, 360], [-1.1, -1.1], color='k')
# 3 JPL picks
plt.plot([176], [-1.2], label='R3 Ia', color=colors[0], marker="D",
         markersize=5)
plt.plot([130, 218], [-1.2, -1.2], color=colors[0])
plt.plot([162], [-1.3], label='R3 Ib', color=colors[0], marker="D",
         markersize=5)
plt.plot([122, 214], [-1.3, -1.3], color=colors[0], linestyle="dashed")
plt.plot([166], [-1.4], label='R3 Ic', color=colors[0], marker="D",
         markersize=5)
plt.plot([120, 220], [-1.4, -1.4], color=colors[0], linestyle="dotted")
# Carrasco
#plt.plot([293], [0.9], label='R3 II', color=colors[1], marker="s",
#         markersize=5)
#plt.plot([266, 301], [0.9, 0.9], color=colors[1])
# DK
plt.plot([136], [-1.5], label='R3 III', color=colors[2], marker='D',
         markersize=5)
plt.plot([7.5, 148.5], [-1.5, -1.5], color=colors[2])
ax.text(270, -1.35, "R3", fontsize=12)
#Add in MQS baz pick
plt.plot([MQSbaz, MQSbaz], [-1.1, -1.6], 'k-.')
ax.add_patch(Rectangle((MQSbaz - 2.*MQSbaz_sigma, -1.6), 4.*MQSbaz_sigma,
                       0.5, color='grey', alpha=0.5))




ax.yaxis.set_ticklabels([])
ax.yaxis.set_ticks([])
ax.xaxis.set_ticks([0, 90, 180, 270, 360])
ax.set_xlim([0, 360])
ax.set_ylim([-1.6, 0.1])
ax.text(labelxoff, 0.96, "E", transform=ax.transAxes, fontsize=16)


# plt.show()
plt.savefig("Figure2_v3.png")
