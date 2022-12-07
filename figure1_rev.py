"""
Make the combined seismic record and picks figure
"""
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
from matplotlib import ticker
from obspy.core import read, UTCDateTime
from obspy.signal.filter import envelope
import numpy as np

# Data file info and reference times
datafile_raw = './data/S1222a.BH.vbbrotate.mseed'
datafile_dg = './data/S1222a.BH.vbbrotate.deglitched.mseed'
origin = UTCDateTime(2022,5,4,23,23,6,570000)
reftime = UTCDateTime(2022,5,4,23,27,45,690000)
comp = 'BHZ'

# Set up the overall grid for all the subplots
fig = plt.figure(figsize=(6,9), dpi=300)
grid = plt.GridSpec(11, 6, wspace=0.3)

axs = [plt.subplot(grid[0, :]), plt.subplot(grid[4:7, :3]),
       plt.subplot(grid[4:7, 3:]), plt.subplot(grid[8:, :3]),
       plt.subplot(grid[8:, 3:]), plt.subplot(grid[2, :2]),
       plt.subplot(grid[2, 2:4]), plt.subplot(grid[2, 4:6])]

labelxoff = [0.01, 0.03, 0.9, 0.03]
labelyoff = [0.1, 0.1, 0.1, 0.9]

# Plot up the reference seismograms
ax = axs[0]
plt.sca(ax)

stz_raw = read(datafile_raw).select(channel=comp)
stz_dg = read(datafile_dg).select(channel=comp)
stime = origin
etime = origin + 9000
f1 = 0.02
f2 = 0.05

stz_raw.filter('bandpass', freqmin=f1, freqmax=f2, corners=2, zerophase=True)
stz_dg.filter('bandpass', freqmin=f1, freqmax=f2, corners=2, zerophase=True)
stz_raw.trim(starttime=stime, endtime=etime)
stz_dg.trim(starttime=stime, endtime=etime)

t_raw = stz_raw[0].times(type="matplotlib")
t_dt_raw = np.array([mdates.num2date(x) for x in t_raw])
t_dg = stz_dg[0].times(type="matplotlib")
t_dt_dg = np.array([mdates.num2date(x) for x in t_dg])

ax.plot(t_dt_raw, stz_raw[0].data, 'r-')
ax.plot(t_dt_dg, stz_dg[0].data, 'k-')
locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
formatter = mdates.ConciseDateFormatter(locator)
yformatter = ticker.ScalarFormatter(useMathText=True)
yformatter.set_scientific(True)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(yformatter)
plt.xlim(t_dt_raw.min(), t_dt_raw.max())
ymin = stz_dg[0].data.min()*1.05
ymax = stz_dg[0].data.max()*1.05
plt.ylim(ymin, ymax)

# Add in some boxes for later plot windows
R1_t1 = UTCDateTime(2022,5,4,23,34,30)
R1_t2 = R1_t1 + 4*60
width = R1_t2.datetime - R1_t1.datetime
ax.add_patch(Rectangle((R1_t1.datetime, ymin), width, ymax-ymin,
                       color='green', alpha=0.5))
R2_t1 = UTCDateTime(2022,5,5,1,12)
R2_t2 = R2_t1 + 4*60
width = R2_t2.datetime - R2_t1.datetime
ax.add_patch(Rectangle((R2_t1.datetime, ymin), width, ymax-ymin,
                       color='green', alpha=0.5))
R3_t1 = UTCDateTime(2022,5,5,1,37)
R3_t2 = R3_t1 + 5*60
width = R3_t2.datetime - R3_t1.datetime
ax.add_patch(Rectangle((R3_t1.datetime, ymin), width, ymax-ymin,
                       color='green', alpha=0.5))
R2_alt_t2 = R2_t1 + 10*60
width = R2_alt_t2.datetime - R2_t2.datetime
ax.add_patch(Rectangle((R2_t2.datetime, ymin), width, ymax-ymin,
                       color='yellow', alpha=0.5))
ax.text(labelxoff[0], labelyoff[0], 'A', transform=ax.transAxes, fontsize=14)

# R1 seismogram zoom
ax = axs[5]
plt.sca(ax)

stz = read(datafile_dg).select(channel=comp)
stz.filter('bandpass', freqmin=0.02, freqmax=0.09, corners=2, zerophase=True)

stime = R1_t1 - 2*60
etime = R1_t2 + 3*60
stz.trim(starttime=stime, endtime=etime)

toffset = stime - reftime

t = stz[0].times('relative')

plt.plot(t+toffset, stz[0].data, 'k-')
ymin = stz[0].data.min() * 1.15
ymax = stz[0].data.max() * 1.15
plt.ylim(ymin, ymax)
plt.xlim(t.min()+toffset, t.max()+toffset)

width = R1_t2 - R1_t1
ax.add_patch(Rectangle((R1_t1-reftime, ymin), width, ymax-ymin, color='green',
                       alpha=0.5))
ax.set_title("R1 zoom")
# ax.set_xlabel("seconds relative to MQS P")
ax.text(labelxoff[1], labelyoff[1], 'B', transform=ax.transAxes, fontsize=14)
ax.set_yticks([])
ax.set_yticklabels([])

# R2 seismogram zoom
ax = axs[6]
plt.sca(ax)

stz = read(datafile_dg).select(channel=comp)
stz.filter('bandpass', freqmin=0.02, freqmax=0.04, corners=2, zerophase=True)

stime = R2_t1 - 2*60
# etime = R3_t2 + 6*60
etime = R2_alt_t2 + 3*60
stz.trim(starttime=stime, endtime=etime)

toffset = stime - reftime

t = stz[0].times('relative')

plt.plot(t+toffset, stz[0].data, 'k-')
ymin = stz[0].data.min() * 1.15
ymax = stz[0].data.max() * 1.15
plt.ylim(ymin, ymax)
plt.xlim(t.min()+toffset, t.max()+toffset)

width = R2_t2 - R2_t1
ax.add_patch(Rectangle((R2_t1-reftime, ymin), width, ymax-ymin, color='green',
                       alpha=0.5))
width = R2_alt_t2 - R2_t2
ax.add_patch(Rectangle((R2_t2-reftime, ymin), width, ymax-ymin, color='yellow',
                       alpha=0.5))
ax.set_title("R2 zoom")
ax.set_xlabel("seconds relative to MQS P")
ax.text(labelxoff[1], labelyoff[1], 'C', transform=ax.transAxes, fontsize=14)
ax.set_yticks([])
ax.set_yticklabels([])

# R3 seismogram zoom
ax = axs[7]
plt.sca(ax)

stz = read(datafile_dg).select(channel=comp)
stz.filter('bandpass', freqmin=0.02, freqmax=0.04, corners=2, zerophase=True)

stime = R3_t1 - 2*60
etime = R3_t2 + 3*60
stz.trim(starttime=stime, endtime=etime)

toffset = stime - reftime

t = stz[0].times('relative')

plt.plot(t+toffset, stz[0].data, 'k-')
ymin = stz[0].data.min() * 1.15
ymax = stz[0].data.max() * 1.15
plt.ylim(ymin, ymax)
plt.xlim(t.min()+toffset, t.max()+toffset)

width = R3_t2 - R3_t1
ax.add_patch(Rectangle((R3_t1-reftime, ymin), width, ymax-ymin, color='green',
                       alpha=0.5))
ax.set_title("R3 zoom")
# ax.set_xlabel("seconds relative to MQS P")
ax.text(labelxoff[1], labelyoff[1], 'D', transform=ax.transAxes, fontsize=14)
ax.set_yticks([])
ax.set_yticklabels([])


# Plot R1 filter banks
ax = axs[1]
plt.sca(ax)

cfreqs = [0.015, 0.018, 0.021, 0.025, 0.03, 0.036, 0.043,
          0.05, 0.06, 0.072, 0.086, 0.1]
cfreqs = np.array(cfreqs)
hwidth = 0.3

stime = UTCDateTime(2022,5,4,23,34,30)
refoffset = stime - reftime
etime = stime + 4*60

SWpicks = [-999, -999, 490.25, 490.84, 501.55, 512.42, 521.27, 529.21,
           536.00, 585.77, 586.76, -999]

st = read(datafile_dg)
# st.trim(starttime=stime, endtime=etime)
print(st)

# comps = ['BHZ', 'BHN', 'BHE']
# comps = ['BHZ']

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

tr = st.select(channel=comp)[0].filter('bandpass', freqmin=cfreqs.min(),
                                       freqmax=cfreqs.max())
toffset = 3*tr.data.max()/len(cfreqs)
pickheight = toffset*3
# toffset = 0.2
# fig = plt.figure()
for i, cfreq in enumerate(cfreqs):
    tr = st.copy().select(channel=comp)[0]
    tr.filter('bandpass', freqmin=cfreq*(1-hwidth),
              freqmax=cfreq*(1+hwidth), corners=2, zerophase=True)
    tr.trim(starttime=stime, endtime=etime)
    ic = i%len(colors)
    plt.plot(tr.times() + refoffset, envelope(tr.data) + i*toffset,
             color=colors[ic])
    # plt.plot(tr.times() + refoffset, tr.data + i*toffset, color=colors[ic],
    #          linestyle='dashed')
    # plt.plot([SWpicks[i], SWpicks[i]],
    #          [i*toffset-pickheight, i*toffset+pickheight],
    #          color=colors[ic], linewidth=3)
    plt.plot([SWpicks[i], SWpicks[i]],
             [i*toffset-pickheight, i*toffset+pickheight],
             color=colors[ic], linewidth=3)

# ax = plt.gca()
ax.set_yticks(np.arange(0, cfreqs.size*toffset, toffset))
ax.set_yticklabels(cfreqs)
ax.set_xlabel('seconds relative to MQS P')
ax.set_ylabel('center frequency (Hz)')
ax.set_xlim(stime-reftime, etime-reftime)
ax.set_ylim(-0.5*pickheight, cfreqs.size*toffset+0.5*pickheight)
ax.set_title("R1 window")
ax.text(labelxoff[3], labelyoff[3], 'E', transform=ax.transAxes, fontsize=14)


# Plot R2 filter banks
ax = axs[2]
plt.sca(ax)

stime = UTCDateTime(2022,5,5,1,12)
refoffset = stime - reftime
etime = stime + 4*60

SWpicks = [-999, -999, 6376.38, 6374.18, 6375.23, 6376.91, 6378.28, 6383.71,
           -999, -999, -999, -999]

st = read(datafile_dg)

tr = st.select(channel=comp)[0].filter('bandpass', freqmin=cfreqs.min(),
                                       freqmax=cfreqs.max())
toffset = 0.75*tr.data.max()/len(cfreqs)
pickheight = toffset*3
# toffset = 0.2
# fig = plt.figure()
for i, cfreq in enumerate(cfreqs):
    tr = st.copy().select(channel=comp)[0]
    tr.filter('bandpass', freqmin=cfreq*(1-hwidth),
              freqmax=cfreq*(1+hwidth), corners=2, zerophase=True)
    tr.trim(starttime=stime, endtime=etime)
    ic = i%len(colors)
    plt.plot(tr.times() + refoffset, envelope(tr.data) + i*toffset,
             color=colors[ic])
    # plt.plot(tr.times() + refoffset, tr.data + i*toffset, color=colors[ic],
    #          linestyle='dashed')
    plt.plot([SWpicks[i], SWpicks[i]],
             [i*toffset-pickheight, i*toffset+pickheight],
             color=colors[ic], linewidth=3)

# ax = plt.gca()
ax.set_yticks(np.arange(0, cfreqs.size*toffset, toffset))
ax.set_yticklabels(cfreqs)
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
ax.set_xlabel('seconds relative to MQS P')
# ax.set_ylabel('center frequency (Hz)')
ax.set_xlim(stime-reftime, etime-reftime)
ax.set_ylim(-0.5*pickheight, cfreqs.size*toffset+0.5*pickheight)
ax.set_title("R2 window")
ax.text(labelxoff[3], labelyoff[3], 'F', transform=ax.transAxes, fontsize=14)

# Plot R3 filter banks
ax = axs[3]
plt.sca(ax)

stime = UTCDateTime(2022,5,5,1,37)
refoffset = stime - reftime
etime = stime + 5*60

SWpicks = [-999, -999, 7872.31, 7863.56, 7864.21, 7882.36, 7882.72, -999,
           -999, -999, -999, -999]

st = read(datafile_dg)

tr = st.select(channel=comp)[0].filter('bandpass', freqmin=cfreqs.min(),
                                       freqmax=cfreqs.max())
toffset = 1*tr.data.max()/len(cfreqs)
pickheight = toffset*3
# toffset = 0.2
# fig = plt.figure()
for i, cfreq in enumerate(cfreqs):
    tr = st.copy().select(channel=comp)[0]
    tr.filter('bandpass', freqmin=cfreq*(1-hwidth),
              freqmax=cfreq*(1+hwidth), corners=2, zerophase=True)
    tr.trim(starttime=stime, endtime=etime)
    ic = i%len(colors)
    plt.plot(tr.times() + refoffset, envelope(tr.data) + i*toffset,
             color=colors[ic])
    # plt.plot(tr.times() + refoffset, tr.data + i*toffset, color=colors[ic],
    #          linestyle='dashed')
    plt.plot([SWpicks[i], SWpicks[i]],
             [i*toffset-pickheight, i*toffset+pickheight],
             color=colors[ic], linewidth=3)

# ax = plt.gca()
ax.set_yticks(np.arange(0, cfreqs.size*toffset, toffset))
ax.set_yticklabels(cfreqs)
ax.set_xlabel('seconds relative to MQS P')
ax.set_ylabel('center frequency (Hz)')
ax.set_xlim(stime-reftime, etime-reftime)
ax.set_ylim(-0.5*pickheight, cfreqs.size*toffset+0.5*pickheight)
ax.set_title("R3 window")
ax.text(labelxoff[3], labelyoff[3], 'G', transform=ax.transAxes, fontsize=14)

# Plot alternate R2 window
ax = axs[4]
plt.sca(ax)

stime = UTCDateTime(2022,5,5,1,12)
refoffset = stime - reftime
etime = stime + 10*60

SWpicks = [-999, -999, 6376.38, 6374.18, 6375.23, 6376.91, 6378.28, 6383.71,
           -999, -999, -999, -999]
SWpicks_late = [-999, -999, 6661.96, 6653.51, 6688.06, 6697.87, 6708.01,
                6726.36, -999, -999, -999, -999]

st = read(datafile_dg)

tr = st.select(channel=comp)[0].filter('bandpass', freqmin=cfreqs.min(),
                                       freqmax=cfreqs.max())
toffset = 0.75*tr.data.max()/len(cfreqs)
pickheight = toffset*3
# toffset = 0.2
# fig = plt.figure()
for i, cfreq in enumerate(cfreqs):
    tr = st.copy().select(channel=comp)[0]
    tr.filter('bandpass', freqmin=cfreq*(1-hwidth),
              freqmax=cfreq*(1+hwidth), corners=2, zerophase=True)
    tr.trim(starttime=stime, endtime=etime)
    ic = i%len(colors)
    plt.plot(tr.times() + refoffset, envelope(tr.data) + i*toffset,
             color=colors[ic])
    # plt.plot(tr.times() + refoffset, tr.data + i*toffset, color=colors[ic],
    #          linestyle='dashed')
    plt.plot([SWpicks[i], SWpicks[i]],
             [i*toffset-pickheight, i*toffset+pickheight],
             color=colors[ic], linewidth=3)
    plt.plot([SWpicks_late[i], SWpicks_late[i]],
             [i*toffset-pickheight, i*toffset+pickheight],
             color=colors[ic], linewidth=1.5, linestyle='dashed')

# ax = plt.gca()
ax.set_yticks(np.arange(0, cfreqs.size*toffset, toffset))
ax.set_yticklabels(cfreqs)
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
ax.set_xlabel('seconds relative to MQS P')
# ax.set_ylabel('center frequency (Hz)')
ax.set_xlim(stime-reftime, etime-reftime)
ax.set_ylim(-0.5*pickheight, cfreqs.size*toffset+0.5*pickheight)

ax.set_title("Alternate R2")
ax.text(labelxoff[3], labelyoff[3], 'H', transform=ax.transAxes, fontsize=14)

# fig.tight_layout()
# plt.show()
plt.savefig("Figure1_rev.png")

