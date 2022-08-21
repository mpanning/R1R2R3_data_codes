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
grid = plt.GridSpec(9, 2, hspace=0, wspace=0.2)
axs = [plt.subplot(grid[:3, 0]), plt.subplot(grid[:3, 1]),
       plt.subplot(grid[4:7, 0]), plt.subplot(grid[4:7, 1]),
       plt.subplot(grid[8,:])]
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
axs[4].legend(h,l, borderaxespad=0, loc=10, ncol=4)
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

# plt.show()
plt.savefig("Figure2_v2.png")
