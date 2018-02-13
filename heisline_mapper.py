import numpy as np
# import os
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.interpolate import sproot, splrep, splev
# from scipy.interpolate import interp1d
# from mpl_toolkits.mplot3d import Axes3D
import pyregion
# from matplotlib import axes
# from matplotlib import cm


class Error(Exception):
    """Base class for other exceptions"""
    pass


class NotAnOption(Error):
    """Raised when the input value is not an option"""
    pass

gridfolder = 'griddata/'
haimagename = 'HaFullSci.fit'
haerrorname = 'HaFullErr.fit'
siiimagename = 'SIIFullSci.fit'
siierrorname = 'SIIFullErr.fit'
contourfileds9 = 'photometry.reg'


# Function that will Read the images
def reader(imname):
    imagehdu = fits.open(imname)
    imagedco = imagehdu[0].data
    toreturn = []
    for iter1 in range(0, len(imagedco)):
        dummy = imagedco[iter1]
        for jter1 in range(0, len(dummy)):
            toreturn.append(dummy[jter1])
    imagehdu.close()
    return toreturn

densityaccept = [0.01, 0.1, 1, 10, 100, 1000]
Baccept = [0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 1, 10]

# Photometric Data etc.
airmass = 1.0
# atmospheric extinction coefficient
k = 0.099
# zero-points
haZP = 4.957
siiZP = 5.077
# signal to noise thresholds
ratiosnthresh = 0.0
hasnthresh = 0.5
siisnthresh = 0.0
# medium density in cm^-3
density = 100
# magnetic field in uG
Bkept = 5

# read data
while True:
    entered = raw_input("Airmass (%.1f):\t" % airmass)
    try:
        airmass = float(entered)
        break
    except ValueError:
        print("Wrong Input")

while True:
    entered = raw_input("Extinction Coefficient (k) (%.3f):\t" % k)
    try:
        k = float(entered)
        break
    except ValueError:
        print("Wrong Input")

while True:
    entered = raw_input("Ha zero-point (%.3f):\t" % haZP)
    try:
        haZP = float(entered)
        break
    except ValueError:
        print("Wrong Input")

while True:
    entered = raw_input("SII zero-point (%.3f):\t" % siiZP)
    try:
        siiZP = float(entered)
        break
    except ValueError:
        print("Wrong Input")

while True:
    entered = raw_input("S/N cutoff for ratio (%.1f):\t" % ratiosnthresh)
    try:
        ratiosnthresh = float(entered)
        break
    except ValueError:
        print("Wrong Input")

while True:
    entered = raw_input("S/N cutoff for Ha (%.1f):\t" % hasnthresh)
    try:
        hasnthresh = float(entered)
        break
    except ValueError:
        print("Wrong Input")

while True:
    entered = raw_input("S/N cutoff for SII (%.1f):\t" % siisnthresh)
    try:
        siisnthresh = float(entered)
        break
    except ValueError:
        print("Wrong Input")

while True:
    entered = raw_input("Density in cm^{-3} (%.0f):\t" % density)
    try:
        density = int(entered)
        if density in densityaccept:
            break
        else:
            raise NotAnOption("Not in list")
    except ValueError:
        print("Wrong Input")
    except NotAnOption:
        print("Density must be one of: ", densityaccept)

while True:
    entered = raw_input("Magnetic Field in uG (%.0f):\t" % Bkept)
    try:
        Bkept = int(entered)
        if Bkept in Baccept:
            break
        else:
            raise NotAnOption("Not in list")
    except ValueError:
        print("Wrong Input")
    except NotAnOption:
        print("Magnetic field must be one of: ", Baccept)


# In[38]:


# Function that read images as 2D
def read2d(imname):
    imhdu = fits.open(imname)
    readdata = imhdu[0].data
    readhead = imhdu[0].header
    imhdu.close()
    toreturn = [readdata, readhead]
    return toreturn


# In[39]:


# Function that will remove zeros from lists
def zerobad(image, error):
    newimage = np.where(image > 0.0, image, 0.0)
    mask = np.divide(newimage, newimage)
    newerror = np.multiply(mask, error)
    toret = [newimage, newerror]
    return toret


# In[40]:


# convert images from counts to flux
def toflux(imlist, kabs, amass, zerop):
    return np.multiply(1, np.multiply((np.power(10.0, (kabs*amass/2.5))),
                                      (np.multiply((np.power(10.0, (zerop/2.5))), imlist))))


# In[41]:


# filter images based on s/n ratio
def snfilter(image1, errorimage1, sn1, image2='none', errorimage2='none', sn2=0.0):
    if image2 == 'none' or errorimage2 == 'none':
        significance = np.divide(image1, errorimage1)
        newimage = np.where(significance > sn1, significance, 0.0)
        newimage = np.divide(newimage, newimage)
        newimage = np.multiply(newimage, image1)
        newerror = np.multiply(newimage, errorimage1)
        toret = [newimage, newerror]
    else:
        significance1 = np.divide(image1, errorimage1)
        mask1 = np.where(significance1 > sn1, significance1, 0.0)
        mask1 = np.divide(mask1, mask1)
        significance2 = np.divide(image2, errorimage2)
        mask2 = np.where(significance2 > sn2, significance2, 0.0)
        mask2 = np.divide(mask2, mask2)
        mask = np.multiply(mask1, mask2)
        newimage1 = np.multiply(mask, image1)
        newerror1 = np.multiply(mask, errorimage1)
        newimage2 = np.multiply(mask, image2)
        newerror2 = np.multiply(mask, errorimage2)
        toret = [newimage1, newerror1, newimage2, newerror2]
    return toret


# In[42]:


# Function that will remove zeros from lists
def maxvaluefilt(image, error, value):
    newimage = np.where(image < value, image, 0.0)
    mask = np.divide(newimage, newimage)
    newerror = np.multiply(mask, error)
    toret = [newimage, newerror]
    return toret


# In[43]:


def minvaluefilt(image, value):
    newimage = np.where(image > value, image, 0.0)
    mask = np.divide(newimage, newimage)
    newimage = np.multiply(mask, newimage)
    return newimage

# read ha science image
hasci = read2d(haimagename)[0]
# read ha error image
haerr = read2d(haerrorname)[0]
# read sii science image
siisci = read2d(siiimagename)[0]
# read sii error image
siierr = read2d(siierrorname)[0]
# convert ha science image to physial units
hasci = toflux(hasci, k, airmass, haZP)
# convert ha error image to physial units
haerr = toflux(haerr, k, airmass, haZP)
# convert sii science image to physial units
siisci = toflux(siisci, k, airmass, siiZP)
# convert sii error image to physial units
siierr = toflux(siierr, k, airmass, siiZP)
# remove bad pixels from ha image
ha3d = zerobad(hasci, haerr)
hasci = ha3d[0]
haerr = ha3d[1]
# remove bad pixels from sii image
sii3d = zerobad(siisci, siierr)
siisci = sii3d[0]
siierr = sii3d[1]
# create pre-S/N-filter ratio image
preratioval = np.divide(siisci, hasci)
# create pre-S/N-filter ratio error image
preratioerr = np.sqrt((np.square(np.divide(siierr, hasci)))+(np.square(np.divide((np.multiply(siisci, haerr)),
                                                                                 (np.square(hasci))))))
# filter pre-S/N-filter ratio image and error for bad values
preratio3d = maxvaluefilt(preratioval, preratioerr, 0.78)
preratioval = preratio3d[0]
preratioerr = preratio3d[1]
# filter ha and sii images based on S/N
snfilt3d = snfilter(image1=hasci, image2=siisci, errorimage1=haerr, errorimage2=siierr, sn1=hasnthresh, sn2=siisnthresh)
hasci = snfilt3d[0]
haerr = snfilt3d[1]
siisci = snfilt3d[2]
siierr = snfilt3d[3]
# create post-S/N-filter ratio image
ratioval = np.divide(siisci, hasci)
# create post-S/N-filter ratio error image
ratioerr = np.sqrt((np.square(np.divide(siierr, hasci)))+(np.square(np.divide((np.multiply(siisci, haerr)),
                                                                              (np.square(hasci))))))
# filter ratio image and error image based on S/N
ratio3d = snfilter(image1=ratioval, errorimage1=ratioerr, sn1=ratiosnthresh)
ratioval = ratio3d[0]
ratioerr = ratio3d[1]
# filter ratio image and error image for bad values
ratio3d = maxvaluefilt(ratioval, ratioerr, 0.78)
ratioval = ratio3d[0]
ratioerr = ratio3d[1]
# heatmap ratio pre- and post- S/N filtering
# '''
# %matplotlib notebook
fig = plt.figure()
preplot = fig.add_subplot(1, 2, 1)
preplot.set_title("[SII]/Ha before S/N filtering")
preplot.imshow(preratioval, interpolation='nearest')
postplot = fig.add_subplot(1, 2, 2)
postplot.set_title("[SII]/Ha after S/N filtering")
postplot.imshow(ratioval, interpolation='nearest')
plt.show()
# '''

# heatmap ratio
plt.clf()
colorbarred = plt.imshow(ratioval, interpolation='nearest')
plt.colorbar()
plt.show()


# In[46]:


# interpolate the SII/Ha vs Shock Velocity curve
# name of the file containing the grid
gridfile = gridfolder+'oiiioverhb_siioverhanii_solar_dens'+str(density)+'_grid1.txt'
# open the file
filer = open(gridfile, 'r')
# write every line to the list
gridlineslist = []
for line in filer:
    gridlineslist.append(line)
# delete the first 11 not-eeded lines
del gridlineslist[:11]
# remove the \n from the list
gridlineslistold = gridlineslist
gridlineslist = []
for i in gridlineslistold:
    gridlineslist.append(i.strip('\n'))
# Read parameters in lists
modelB = []
modelvel = []
modelSII = []
modelNIInHa = []
modelHa = []
for item1 in gridlineslist:
    dummy1 = item1.split(" ")
    dummy2 = []
    for item2 in dummy1:
        if item2 != '':
            dummy2.append(float(item2))
    modelB.append(dummy2[0])
    modelvel.append(dummy2[1])
    modelSII.append(dummy2[2])
    modelNIInHa.append(dummy2[3])
    modelHa.append(dummy2[5])
modelSIIoverHa = np.divide(modelSII, modelHa)
modelSIIoverNIInHa = np.divide(modelSII, modelNIInHa)
# Isolate the needed magnetic field
modelratio = []
modelvelocity = []
for i in range(0, len(modelB)):
    if modelB[i] == Bkept:
        modelratio.append(modelSIIoverNIInHa[i])
        modelvelocity.append(modelvel[i])
plt.clf()
plt.plot(modelvelocity, modelratio, '.')
modelvelocityinterp = np.linspace(min(modelvelocity), max(modelvelocity), num=1000, endpoint=True)
modelratiotck = splrep(modelvelocity, modelratio)
modelratiointerp = splev(modelvelocityinterp, modelratiotck)
plt.plot(modelvelocityinterp, modelratiointerp)
plt.title(r"$\frac{[SII]}{H\alpha+[NII]}$ Diagnostic", y=1.02)
plt.xlabel(r"Velocity $(\frac{km}{s})$")
plt.ylabel(r"$\frac{[SII]}{H\alpha+[NII]}$")
plt.show()

rationew = ratioval[np.logical_not(np.isnan(ratioval))]
rationew = minvaluefilt(rationew, 0.0)
rationew = rationew[np.logical_not(np.isnan(rationew))]
uppercap = np.max(modelratio)
lowercap = np.min(modelratio)
print("Max model ratio: %.6f" % uppercap)
print("Min model ratio: %.6f" % lowercap)

np.set_printoptions(threshold=np.nan)
# print(ratioval)
# print(ratioval[20])
# print(ratioval[20,6])

ratioplus = np.add(ratioval, ratioerr)
ratiominus = np.subtract(ratioval, ratioerr)
# ratiominus=np.where(ratiominus>=0.0, ratiominus, np.nan)
ratioval = np.where(ratioval >= 0.000000000001, ratioval, np.nan)
plt.clf()
colorbarred = plt.imshow(ratioval, interpolation='nearest')
plt.colorbar()
plt.show()

plt.clf()
colorbarred = plt.imshow(ratioplus, interpolation='nearest')
plt.colorbar()
plt.show()

plt.clf()
colorbarred = plt.imshow(ratiominus, interpolation='nearest')
plt.colorbar()
plt.show()

dimratioval = np.shape(ratioval)
dimratioplus = np.shape(ratioplus)
dimratiominus = np.shape(ratiominus)
velocity = np.empty([2, dimratioval[0], dimratioval[1]])
velocityplus = np.empty([2, dimratioplus[0], dimratioplus[1]])
velocityminus = np.empty([2, dimratiominus[0], dimratiominus[1]])


def determinevelocity(ratiotosub):
    if ratiotosub >= uppercap:
        ratiotosub = uppercap
    if ratiotosub <= lowercap:
        ratiotosub = lowercap
    if np.isnan(ratiotosub):
        return np.nan
    newmodelratiotck = splrep(modelvelocity, np.subtract(modelratio, ratiotosub))
    velocities = sproot(newmodelratiotck)
    velocitiesret = [np.min(velocities), np.max(velocities)]
    return velocitiesret


def determinevelocitycap(ratiotosub):
    if ratiotosub >= uppercap:
        return [0.0, 0.0]
    if ratiotosub <= lowercap:
        return [0.0, 0.0]
    if np.isnan(ratiotosub):
        return np.nan
    newmodelratiotck = splrep(modelvelocity, np.subtract(modelratio, ratiotosub))
    velocities = sproot(newmodelratiotck)
    velocitiesret = [np.min(velocities), np.max(velocities)]
    return velocitiesret

for i in range(0, dimratioval[0]):
    for j in range(0, dimratioval[1]):
        pixelvel = determinevelocity(ratioval[i, j])
        if np.any(np.isnan(pixelvel)):
            velocity[0, i, j] = np.nan
            velocity[1, i, j] = np.nan
        else:
            velocity[0, i, j] = pixelvel[0]
            velocity[1, i, j] = pixelvel[1]
for i in range(0, dimratioplus[0]):
    for j in range(0, dimratioplus[1]):
        pixelvel = determinevelocitycap(ratioplus[i, j])
        if np.any(np.isnan(pixelvel)):
            velocityplus[0, i, j] = np.nan
            velocityplus[1, i, j] = np.nan
        else:
            velocityplus[0, i, j] = pixelvel[0]
            velocityplus[1, i, j] = pixelvel[1]
for i in range(0, dimratiominus[0]):
    for j in range(0, dimratiominus[1]):
        pixelvel = determinevelocitycap(ratiominus[i, j])
        if np.any(np.isnan(pixelvel)):
            velocityminus[0, i, j] = np.nan
            velocityminus[1, i, j] = np.nan
        else:
            velocityminus[0, i, j] = pixelvel[0]
            velocityminus[1, i, j] = pixelvel[1]

plt.clf()
f_xray = fits.open(haimagename)
contourloaded = pyregion.open(contourfileds9).as_imagecoord(f_xray[0].header)
patch_list, text_list = contourloaded.get_mpl_patches_texts()
colorbarred = plt.subplot(111)
tobar = colorbarred.imshow(velocity[0], interpolation='nearest', origin='lower')
for p in patch_list:
    colorbarred.add_patch(p)
for t in text_list:
    colorbarred.add_artist(t)
# plt.colorbar(ax=tobar)
plt.title("Velocity 1")
plt.xlabel("XPixel")
plt.ylabel("YPixel")
plt.show()

plt.clf()
colorbarred = plt.imshow(velocity[0], interpolation='nearest')
plt.colorbar()
plt.title("Velocity 1")
plt.xlabel("XPixel")
plt.ylabel("YPixel")
plt.show()

plt.clf()
colorbarred = plt.imshow(velocityplus[0], interpolation='nearest')
plt.colorbar()
plt.title("Velocity 1 Upper Bound @ 1sigma")
plt.xlabel("XPixel")
plt.ylabel("YPixel")
plt.show()

plt.clf()
colorbarred = plt.imshow(velocityminus[0], interpolation='nearest')
plt.colorbar()
plt.title("Velocity 1 Lower Bound @ 1sigma")
plt.xlabel("XPixel")
plt.ylabel("YPixel")
plt.show()

plt.clf()
f_xray = fits.open(haimagename)
contourloaded = pyregion.open(contourfileds9).as_imagecoord(f_xray[0].header)
patch_list, text_list = contourloaded.get_mpl_patches_texts()
colorbarred = plt.subplot(111)
tobar = colorbarred.imshow(velocity[1], interpolation='nearest', origin='lower')
for p in patch_list:
    colorbarred.add_patch(p)
for t in text_list:
    colorbarred.add_artist(t)
# plt.colorbar(ax=tobar)
plt.title("Velocity 2")
plt.xlabel("XPixel")
plt.ylabel("YPixel")
plt.show()

plt.clf()
colorbarred = plt.imshow(velocity[1], interpolation='nearest')
plt.colorbar()
plt.title("Velocity 2")
plt.xlabel("XPixel")
plt.ylabel("YPixel")
plt.show()

plt.clf()
colorbarred = plt.imshow(velocityplus[1], interpolation='nearest')
plt.colorbar()
plt.title("Velocity 2 Upper Bound @ 1sigma")
plt.xlabel("XPixel")
plt.ylabel("YPixel")
plt.show()

plt.clf()
colorbarred = plt.imshow(velocityminus[1], interpolation='nearest')
plt.colorbar()
plt.title("Velocity 2 Lower Bound @ 1sigma")
plt.xlabel("XPixel")
plt.ylabel("YPixel")
plt.show()
