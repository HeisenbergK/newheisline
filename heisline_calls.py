from pyraf import iraf
import sys
import os
import time
import numpy as np
from astropy.io import ascii, fits
import scipy.spatial
import pyregion

heislineversion = 3.1
date = "February 19 2018"


# check for interactivity
def masterinteractivitycheck(typed):
    messages = ["Do you want interactive mode (to be able to replace files on fly)? (Y/n) \t",
                "If you want replace BIAS.fit do it now then if you want to turn off interactivity answer n \t",
                "If you want replace the bias-subtracted images (*_b.fit) do it now then if you want to turn "
                "off interactivity answer n \t",
                "If you want replace the flat images (FLAT_*_norm.fit) do it now then if you want to turn "
                "off interactivity answer n \t",
                "If you want replace the flat-fielded images (*_f.fit) do it now then "
                "if you want to turn off interactivity answer n \t",
                "If you want replace the reference image for the alignment do it now editing the merged_b_f and "
                "merged_b_f_sh lists putting the desired image first\t",
                "If you want replace the coordinates file (aligcoordslog)\t"]
    response = 'y'
    while True:
        response = raw_input(messages[typed])
        if response == 'n' or response == 'N':
            response = 'n'
            break
        elif response == 'y' or response == 'Y' or response == '':
            response = 'y'
            break
    return response


# ask whether user want something
def ask(typed):
    messages = ["Do you want bias-subtraction?\t",
                "Do you want flat-fielding?\t",
                "Do you want to align images?\t",
                "Do you want to median combine the images?\t",
                "Do you want to add WCS to the images?\t",
                "Do you want to calculate the scaling factor of the continuum images?\t",
                "Do you want to scale-subtract the continuum images from the narrow-band ones?\t",
                "Do you want to follow the automated method to flux calibrate (y) or input the extinction "
                "coefficient and Zero Point yourself (n)?\t",
                "Do you want to pivot the atmospheric extinction coefficient?\t",
                "Do you want to fluxcalibrate using IRAF (only if you now photutils will fail)?\t"]
    response = 'y'
    while True:
        response = raw_input(messages[typed])
        if response == 'n' or response == 'N':
            response = 'n'
            break
        if response == 'y' or response == 'Y':
            response = 'y'
            break
    return response


# average-combines bias frames to master bias frame
def combinebias(biaslist):
    iraf.unlearn('zerocombine')
    iraf.noao.imred.ccdred.zerocombine(biaslist, output="BIAS.fit", combine="average", rdnoise=8.1, gain=2.867,
                                       ccdtype="")
    iraf.unlearn('zerocombine')


# user chooses filters to be used
def choosefilters():
    numfil = 0
    while True:
        numfil = input("Number of filters to be used (counting continuum):\t")
        if isinstance(numfil, int):
            break
    numfil = int(numfil)
    filtnames = []
    i = 0
    while True:
        if i < numfil:
            tmp = raw_input("Enter name of filter %01d:\t" % (i + 1))
            filtnames.append(tmp)
            i += 1
        else:
            break
    print(filtnames)
    return filtnames


# subtracts bias from files
def subtractbias(imlist, biasimage):
    iraf.unlearn("ccdproc")
    iraf.noao.imred.ccdred.ccdproc(images=("@" + imlist), output=("@" + imlist + "_b"), ccdtype="", fixpix="no",
                                   overscan="no", trim="no", darkcor="no", flatcor="no", zero=biasimage)
    iraf.unlearn("ccdproc")


# median-combines flat frames to master flat frame for filter
def combineflats(flatlist, filtname):
    iraf.unlearn("immatch")
    iraf.images.immatch.imcombine(input=("@" + flatlist + "_b"), output=("FLAT_" + filtname + ".fit"),
                                  combine="median", rdnoise=8.1, gain=2.867)
    iraf.unlearn("immatch")
    iraf.unlearn("imstatistics")
    flat_mean = iraf.images.imutil.imstatistics(images=("FLAT_" + filtname + ".fit"), fields="mean", format="No",
                                                Stdout=1)
    iraf.unlearn("imstatistics")
    flat_mean_num = float(flat_mean[0])
    if flat_mean_num == 0:
        print("Error on flat image %s" % flatlist)
        sys.exit()
    iraf.unlearn("imarith")
    iraf.images.imutil.imarith(operand1=("FLAT_" + filtname + ".fit"), op="/", operand2=flat_mean_num,
                               result=("FLAT_" + filtname + "_norm.fit"))
    iraf.unlearn("imarith")


# divides images by normalized flat frame
def flatfielding(imagelist, filtname):
    iraf.unlearn("imarith")
    iraf.images.imutil.imarith(operand1=("@" + imagelist + "_b"), op="/", operand2=("FLAT_" + filtname + "_norm.fit"),
                               result=("@" + imagelist + "_b_f"))
    iraf.unlearn("imarith")


# aligns images using imalign
def aligning(referenceim):
    print("Now the images will be shifted")
    bbox = input("Enter a rough estimate of the max shift needed\t")
    sbox = input("Enter a rough estimate of the FWHM\t")
    bbox = int(2 * bbox)
    sbox += 2
    os.system("ds9 &")
    print("You are now to select 5-10 stars that are isolated and well-defined")
    referenceim = referenceim.rstrip("\n")
    time.sleep(10)
    iraf.unlearn("display")
    iraf.images.tv.display(image=referenceim, frame=1, fill="Yes")
    iraf.unlearn("display")
    iraf.unlearn("imexamine")
    iraf.images.tv.imexamine(input=referenceim, frame=1, logfile="aligcoordslog", keeplog="Yes")
    iraf.unlearn("imexamine")
    masterinteractivitycheck(6)
    correct = "n"
    counter = 1
    while True:
        if correct == "n":
            sbox -= 1
            bbox *= 1.1
            if counter > 1:
                os.system("rm -rf *_sh.fit")
            iraf.unlearn("imalign")
            iraf.images.immatch.imalign(input="@merged_b_f", reference=referenceim, coords="aligcoordslog",
                                        output="@merged_b_f_sh", niterate=200, boxsize=sbox, bigbox=bbox)
            iraf.unlearn("imalign")
        print("Please examine if the alignment was correct. Close ds9 and answer y/n. If you answer with something "
              "other than 'y' or 'n' the result won't change")
        os.system("ds9 *_sh.fit")
        correct = raw_input()
        if counter == 10:
            break
        if correct == "y":
            break
        counter += 1
    if counter == 10:
        failalign = raw_input("I cannot align the images. Either the max shift and FWHM was incorrect or the images "
                              "are bad.  Please align them yourself, replace the existing *_sh.fit files aand "
                              "hit enter")
        print(failalign)


def combining(imagelist, filtname, basename):
    iraf.unlearn("imcombine")
    iraf.images.immatch.imcombine(input=("@" + imagelist), output=(basename + "_" + filtname + ".fit"),
                                  combine="median", rdnoise=8.1, gain=2.867)
    iraf.unlearn("imcombine")


def astrometry(directory, filtname, basename, roughestimatera, roughestimatedec):
    os.chdir(directory + "/" + filtname)
    os.system("solve-field --ra=%s --dec=%s -5 0.5 --out=%s --sigma 25 -t=2 %s" %
              (roughestimatera, roughestimatedec, "test.fit", basename + "_" + filtname + ".fit"))
    os.system("mv test.new %s" % (basename + "_" + filtname + "_W.fit"))
    os.system("rm %s %s %s %s %s %s %s %s %s %s" %
              ("test-ngc.png", "test.wcs", "test-objs.png", "test-indx.png", "test-indx.xyls", "test.rdls",
               "test.axy", "test.solved", "test.match", "test.corr"))
    if not os.path.isfile(basename + "_" + filtname + "_W.fit"):
        os.system("cp %s %s" % ((basename + "_" + filtname + ".fit"), (basename + "_" + filtname + "_W.fit")))


def scaling(directory, filtnames, basename):
    os.system("mkdir scaling")
    os.chdir("scaling")
    for i in filtnames:
        os.system("cp %s/%s/photometry1.mag ./phot%s" % (directory, i, i))
        os.system("cp %s/%s/%s ." % (directory, i, (basename + "_" + i + "_W_allstarfin.fit")))
    numnarr = 0
    while True:
        print(filtnames)
        numnarr = input("Out of the shown filters, how many are narrowband (not continuum)?\t")
        if isinstance(numnarr, int):
            break
    narrows = []
    i = 0
    while True:
        if i < numnarr:
            tmpnarrname = raw_input("What is the name of the #%01d narrow filter name?\t" % (i + 1))
            narrows.append(tmpnarrname)
            i += 1
        else:
            break
    pairs = []
    i = 0
    while True:
        if i < numnarr:
            tmpcontname = raw_input(
                "What is the name of the continuum filter corresponding to the %s narrow filter? \t" % (narrows[i]))
            tmppair = [narrows[i], tmpcontname]
            pairs.append(tmppair)
            i += 1
        else:
            break
    print(pairs)
    scfact = []
    scfactstd = []
    for i in range(0, len(pairs)):
        currentpair = pairs[i]
        currentnarr = currentpair[0]
        currentcont = currentpair[1]
        narrowcatalogfile = ("phot%s" % currentnarr)
        continuumcatalogfile = ("phot%s" % currentcont)
        narrowcatalog = ascii.read(narrowcatalogfile, format='daophot', include_names=['XCENTER', 'YCENTER', 'FLUX'])
        continuumcatalog = ascii.read(continuumcatalogfile, format='daophot',
                                      include_names=['XCENTER', 'YCENTER', 'FLUX'])
        xcennarr = narrowcatalog['XCENTER']
        ycennarr = narrowcatalog['YCENTER']
        xcencont = continuumcatalog['XCENTER']
        ycencont = continuumcatalog['YCENTER']
        prelimnarr = [xcennarr, ycennarr]
        prelimcont = [xcencont, ycencont]
        prelimnarr = np.transpose(prelimnarr)
        prelimcont = np.transpose(prelimcont)
        narrtree = scipy.spatial.cKDTree(prelimnarr, leafsize=100)
        excludenarr = []
        for item in prelimnarr:
            excludenarr.append(narrtree.query(item, k=1, distance_upper_bound=1.2))

        conttree = scipy.spatial.cKDTree(prelimcont, leafsize=100)
        excludecont = []
        for item in prelimcont:
            excludecont.append(conttree.query(item, k=1, distance_upper_bound=1.2))

        j = len(excludenarr) - 1
        while j >= 0:
            tupleread = excludenarr[j]
            if tupleread[0] == 0.0:
                del excludenarr[j]
            j -= 1

        j = len(excludecont) - 1
        while j >= 0:
            tupleread = excludecont[j]
            if tupleread[0] == 0.0:
                del excludecont[j]
            j -= 1
        contdone = []

        for j in range(0, len(prelimcont)):
            result = narrtree.query(prelimcont[j], k=1, distance_upper_bound=1.2)
            result = result + (j,)
            contdone.append(result)
        j = len(contdone) - 1
        while j >= 0:
            tupleread = contdone[j]
            if np.isinf(tupleread[0]):
                del contdone[j]
            j -= 1
        ratio = []
        for j in range(0, len(contdone)):
            tupleread = contdone[j]
            fluxnarrpoint = narrowcatalog['FLUX'][tupleread[1]]
            fluxcontpoint = continuumcatalog['FLUX'][tupleread[2]]
            ratiocur = fluxnarrpoint / fluxcontpoint
            if (not np.isnan(ratiocur)) and (not np.isinf(ratiocur)):
                ratio.append(ratiocur)
        scfact.append(np.median(ratio))
        scfactstd.append(np.std(ratio))
    packs = []
    for i in range(0, len(pairs)):
        currentpack = pairs[i]
        currentpack.append(scfact[i])
        currentpack.append(scfactstd[i])
        packs.append(currentpack)
    scalefile = open('scalefactors', 'w')
    for item in packs:
        scalefile.write("%s,%s,%.5f,%.7f" % (item[0], item[1], item[2], item[3]))
    scalefile.close()
    print(packs)
    return packs, narrows


def scalesub(packs, basename):
    for item in packs:
        narrfilt = item[0]
        contfilt = item[1]
        factorstring = item[2]
        uncertaintystring = item[3]
        factor = float(factorstring)
        uncertainty = float(uncertaintystring)
        factorplus = factor + uncertainty
        factorminus = factor - uncertainty
        narrimagename = (basename + "_" + narrfilt + "_W_allstarfin.fit")
        contimagename = (basename + "_" + contfilt + "_W_allstarfin.fit")
        scaledcontimagename = (basename + "_" + contfilt + "scaledfor" + narrfilt + "_W_allstarfin.fit")
        underscaledcontimagename = (basename + "_" + contfilt + "underscaledfor" + narrfilt + "_W_allstarfin.fit")
        overscaledcontimagename = (basename + "_" + contfilt + "overscaledfor" + narrfilt + "_W_allstarfin.fit")
        finalimagenamezersig = (basename + "_" + narrfilt + "_W_allstarfin_scalesub0.fit")
        finalimagenameplusig = (basename + "_" + narrfilt + "_W_allstarfin_scalesub+sigma.fit")
        finalimagenameminsig = (basename + "_" + narrfilt + "_W_allstarfin_scalesub-sigma.fit")
        iraf.unlearn("imarith")
        # zero sigma
        iraf.imarith(operand1=contimagename, op="*", operand2=factor, result=scaledcontimagename)
        iraf.imarith(operand1=narrimagename, op="-", operand2=scaledcontimagename, result=finalimagenamezersig)
        # factor + 1 sigma
        iraf.imarith(operand1=contimagename, op="*", operand2=factorplus, result=overscaledcontimagename)
        iraf.imarith(operand1=narrimagename, op="-", operand2=overscaledcontimagename, result=finalimagenameplusig)
        # factor - 1 sigma
        iraf.imarith(operand1=contimagename, op="*", operand2=factorminus, result=underscaledcontimagename)
        iraf.imarith(operand1=narrimagename, op="-", operand2=underscaledcontimagename, result=finalimagenameminsig)


def skyeval(imagename, skyregname):
    image = fits.open(imagename)
    image = image[0]
    region = pyregion.open(skyregname).as_imagecoord(image.header)
    mask = region.get_mask(shape=image.data.shape)
    masked = np.multiply(image.data, mask)
    masked = np.array(masked)
    skyvals = []
    for i in range(0, len(masked)):
        for j in range(0, len(masked[i])):
            if masked[i][j] != 0:
                skyvals.append(masked[i][j])
    sky = np.median(skyvals)
    skyerr = np.std(skyvals)
    return sky, skyerr


def patch_up(imagename, problemregion, outim):
    image = fits.open(imagename)
    image = image[0]
    oldheader = image.header
    region = pyregion.open(problemregion).as_imagecoord(image.header)
    mask = region.get_mask(shape=image.data.shape)
    mask = np.subtract(1, mask)
    masked = np.multiply(image.data, mask)
    oldheader.add_history("User regions patched")
    hdu = fits.PrimaryHDU(masked, header=oldheader)  # create fits HDU from binned image
    hdu.writeto(outim)


def binimage(imname, binsize, binnedname):
    binsizestr = str(binsize)
    fitsfile = fits.open(imname)
    im = np.array(fitsfile[0].data)
    xedge = np.shape(im)[0] % binsize
    yedge = np.shape(im)[1] % binsize
    im = im[int(xedge):, int(yedge):]
    binim = np.reshape(im, (np.shape(im)[0]/binsize, binsize, np.shape(im)[1]/binsize, binsize))
    binim = np.sum(binim, axis=3)
    binim = np.sum(binim, axis=1)
    oldheader = fitsfile[0].header  # get header
    oldheader.add_history("Binned in %s x %s bins" % (binsizestr, binsizestr))
    hdu = fits.PrimaryHDU(binim, header=oldheader)  # create fits HDU from binned image
    hdu.writeto(binnedname)


def instphot(imname, sigma, skyval, skyerr, errimname="None"):
    if errimname == "None":
        errimname = imname
    errim = fits.open(errimname)
    errim = np.array(errim[0].data)
    image = fits.open(imname)
    header = image[0].header
    if "EXPTIME" in header:
        expo = header["EXPTIME"]
        expo = int(expo)
    elif "EXPOSURE" in header:
        expo = header["EXPTIME"]
        expo = int(expo)
    else:
        expo = raw_input("Enter the exposure time of the image %s" % str(imname))
        expo = int(expo)

    if "GAIN" in header:
        epadu = header["GAIN"]
        epadu = float(epadu)
    else:
        epadu = raw_input("Enter the gain of the image %s" % str(imname))
        epadu = float(epadu)

    expo = float(expo)
    print("Exposure time: %.2f s, Gain: %.2f epadu" % (expo, epadu))
    image = np.array(image[0].data)
    if (len(errim) != len(image)) or (len(errim[0]) != len(image[0])):
        errimname = imname
    errim = fits.open(errimname)
    errim = np.array(errim[0].data)
    imcount = []
    errcount = []
    for i in range(0, len(image)):
        for j in range(0, len(image[i])):
            imcount.append(image[i][j])
            errcount.append(errim[i][j])
    sign = np.divide(imcount, np.sqrt(np.abs(errcount)))
    counts = []
    for i in range(0, len(sign)):
        if sign[i] >= sigma:
            counts.append(imcount[i])
    bins = len(counts)
    totcount = sum(counts)
    print("Counts Total= %s BinsTotal= %s Sky= %s" % (str(totcount), str(bins), str(skyval)))
    flux = totcount - (bins * skyval)
    error = np.sqrt((np.divide(flux, epadu)) + (bins * skyerr * skyerr))
    minstr = (2.5 * np.log10(expo)) - (2.5 * np.log10(flux))
    merror = np.sqrt(np.square(1.0857 * (np.divide(error, flux))))
    print("Flux=%s Error=%s" % (str(flux), str(error)))
    return minstr, merror


def narrowtot(narrim, narrskyreg, contim, contskyreg, k, zp, sigma, binsize, f, narrbinned, contbinned, narrfinname,
              narrfinerrname):
    narrskyval, narrskyerr = skyeval(narrim, narrskyreg)
    contskyval, contskyerr = skyeval(contim, contskyreg)

    curdir = os.getcwd()
    print("Now the %s image will be opened. Mark problematic regions on the image and save it as a .reg file in the "
          "folder where we are now (%s). Please use the ds9-suggested coordinates." % (narrim, curdir))
    os.system("ds9 %s" % narrim)
    dummy = raw_input("Enter the name of the file you just saved:\t")
    os.system("mv %s problemsnarr.reg" % dummy)
    print("Now the %s image will be opened. Mark problematic regions on the image and save it as a .reg file in the "
          "folder where we are now (%s). Please use the ds9-suggested coordinates." % (contim, curdir))
    os.system("ds9 %s" % contim)
    dummy = raw_input("Enter the name of the file you just saved:\t")
    os.system("mv %s problemscont.reg" % dummy)

    patch_up(narrim, "problemsnarr.reg", "narrpatched.fit")
    patch_up(contim, "problemscont.reg", "contpatched.fit")
    binimage("narrpatched.fit", binsize, narrbinned)
    binimage("contpatched.fit", binsize, contbinned)
    narrdata = fits.open(narrbinned)
    narrhead = narrdata[0].header
    if "AIRMASS" in narrhead:
        airmass = narrhead["AIRMASS"]
        airmass = float(airmass)
    else:
        airmass = raw_input("Enter the airmass of the image %s" % str(narrim))
        airmass = float(airmass)
    narrdata = np.array(narrdata[0].data)
    contdata = fits.open(contbinned)
    contdata = np.array(contdata[0].data)
    asq = np.power(binsize, 2)
    fasq = f * asq
    fsq = np.power(f, 2)
    facta = np.multiply(asq, narrskyval)
    factb = np.multiply(f, contdata)
    factc = np.multiply(fasq, contskyval)
    factd = np.power((np.multiply(asq, narrskyerr)), 2)
    facte = np.multiply(fsq, contdata)
    factf = np.power((np.multiply(fasq, contskyerr)), 2)
    narrfin = np.add((np.subtract(narrdata, (np.add(facta, factb)))), factc)
    narrfinerr = np.sqrt(np.add(narrdata, np.add(factd, np.add(facte, factf))))
    narrskystring = str(narrskyval)
    factorstring = str(f)
    narrhead.add_history("Sky(=%s) and continuum (f=%s) subtracted" % (narrskystring, factorstring))
    hdu = fits.PrimaryHDU(narrfin, header=narrhead)
    hdu.writeto(narrfinname)
    errhdu = fits.PrimaryHDU(narrfinerr)
    errhdu.writeto(narrfinerrname)
    # binnedskyerr = np.multiply(asq, narrskyerr)
    instrmag, instrerr = instphot(narrfinname, sigma, 0, (asq*narrskyerr), narrfinerrname)
    print (instrmag, instrerr)
    mag = zp - np.multiply(airmass, k) + instrmag
    magerr = instrerr
    return mag, magerr
