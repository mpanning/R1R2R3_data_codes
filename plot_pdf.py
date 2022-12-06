"""
Plot up some MQS pdfs
"""

# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import ticker
from scipy.optimize import curve_fit
from obspy.core import UTCDateTime
import xml.etree.ElementTree as ET
import re

quakemlfile = './data/S1222a_1205.xml'

# Open QuakeML file and read the relevant pdfs
tree = ET.parse(quakemlfile)
root = tree.getroot()

# Find all elements with a pdf child and do stuff
ns = 'http://quakeml.org/xmlns/singlestation/1.0'
variable = []
probability = []
tag = []
for elem in root.findall(".//{"+ns+"}pdf/.."):
    tag.append(re.split('{|}',elem.tag)[-1]) # get rid of namespace
    print(tag[-1])
    variable.append(np.array(elem.findall(".//{"+ns+"}variable")[0].text.split()))
    probability.append(np.array(elem.findall(".//{"+ns+"}probability")[0].text.split()))

tag = np.array(tag)
indx = np.where(tag == 'azimuth')[0][0]
baz = np.array(variable[indx], dtype=np.dtype('float_'))
bazpdf = np.array(probability[indx], dtype=np.dtype('float_'))

def gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

popt, pcov = curve_fit(gauss, baz, bazpdf, p0=[0.0001, 102, 10])

plt.plot(baz, bazpdf)
plt.plot(baz, gauss(baz, *popt))
ax = plt.gca()
ax.text(0.7, 0.9, "Max: {:.0f}".format(baz[np.argmax(bazpdf)]),
        transform=ax.transAxes, fontsize=12)
ax.text(0.7, 0.85, "Mean: {:.2f}".format(popt[1]), transform=ax.transAxes,
        fontsize=12)
ax.text(0.7, 0.8, "$\sigma$: {:.2f}".format(popt[2]), transform=ax.transAxes,
        fontsize=12)
ax.set_xlabel("Backazimuth (degrees east of north)")
# ax.set_ylabel("Probability density")

plt.savefig("MQS_baz.png")
# plt.show()


plt.clf()

indx = np.where(tag == 'distance')[0][0]
dist = np.array(variable[indx], dtype=np.dtype('float_'))
distpdf = np.array(probability[indx], dtype=np.dtype('float_'))

popt, pcov = curve_fit(gauss, dist, distpdf, p0=[150, 35, 1])

plt.plot(dist, distpdf)
plt.plot(dist, gauss(dist, *popt))

ax = plt.gca()
ax.text(0.7, 0.9, "Max: {:.1f}".format(dist[np.argmax(distpdf)]),
        transform=ax.transAxes, fontsize=12)
ax.text(0.7, 0.85, "Mean: {:.2f}".format(popt[1]), transform=ax.transAxes,
        fontsize=12)
ax.text(0.7, 0.8, "$\sigma$: {:.2f}".format(popt[2]), transform=ax.transAxes,
        fontsize=12)
ax.set_xlabel("Distance (degrees)")
ax.set_xlim(0, 90)

# plt.show()
plt.savefig("MQS_dist.png")


plt.clf()

indx = np.where(tag == 'originTime')[0][0]
t0_from_xml = variable[indx]
t0 = np.array([UTCDateTime(x) for x in t0_from_xml])
t0_mdates = np.array([mdates.num2date(x.matplotlib_date) for x in t0])
t0_index = np.arange(len(t0))
t0pdf = np.array(probability[indx], dtype=np.dtype('float_'))

t0_max = t0[np.argmax(t0pdf)]
dt0 = t0[np.argmax(t0pdf)]-t0[np.argmax(t0pdf)-1]

popt, pcov = curve_fit(gauss, t0_index, t0pdf, p0=[0.1, np.argmax(t0pdf), 20])
meantime = t0[int(popt[1])] + dt0*(popt[1]-int(popt[1]))

plt.plot(t0_mdates, t0pdf)
plt.plot(t0_mdates, gauss(t0_index, *popt))
ax = plt.gca()
locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
formatter = mdates.ConciseDateFormatter(locator)
yformatter = ticker.ScalarFormatter(useMathText=True)
yformatter.set_scientific(True)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(yformatter)

ax.text(0.55, 0.9,
        "Max: {}".format(t0[np.argmax(t0pdf)].strftime("%H:%M:%S.%f")),
        transform=ax.transAxes, fontsize=12)
ax.text(0.55, 0.85, "Mean: {}".format(meantime.strftime("%H:%M:%S.%f")),
        transform=ax.transAxes, fontsize=12)
ax.text(0.55, 0.8, "$\sigma$: {:.2f}".format(popt[2]*dt0),
        transform=ax.transAxes, fontsize=12)
ax.set_xlabel("t$_0$ (UTC)")
# ax.set_xlim(0, 90)

# plt.show()
plt.savefig("MQS_t0.png")


             
