import os
import numpy as np
import time
import progressbar
from heisline_readers import *

heislineversion = 2.8
date = "February 13 2018"

def subtractskyandbin(ralln, rpren, haalln, hapren, siialln, siipren, problemfilen, contourfilen, skyvalfilen, rharatio, rsiiratio, binsize, imdir):
    os.chdir(imdir)
    allimages = [ralln, haalln, siialln]
    preimages = [rpren, hapren, siipren]

    print("Calculating R sky:" )
    rtosky = 'byprod/'+ralln.strip('.fit')+'_tosky.fit'
    os.system('dmcopy "'+ralln+'[(#1,#2)=region(' + skyvalfilen + ')]" '+rtosky)
    rtoskylist = reader(rtosky)
    rtoskylistmod = removezeros(rtoskylist)
    rskyval = np.median(rtoskylistmod)
    rskyerr = np.std(rtoskylistmod)
    print("DONE!")

    rallsksbn = 'byprod/'+ralln.strip('.fit')+'_sksub.fit'
    haallrsub = 'byprod/'+haalln.strip('.fit')+'_rsub.fit'
    siiallrsub = 'byprod/'+siialln.strip('.fit')+'_rsub.fit'
    haallrsubmaskedn = 'byprod/'+haalln.strip('.fit')+'_rsubmasked.fit'
    siiallrsubmaskedn = 'byprod/'+siialln.strip('.fit')+'_rsubmasked.fit'
    readr = read2d(ralln)
    rdata = readr[0]
    rhead = readr[1]
    rdatanew = np.subtract(rdata, rskyval)
    fits.writeto(rallsksbn, rdatanew, header=rhead)
    rdata = []
    rhead = []
    readr = []
    rdatanew = []

    readr = read2d(rallsksbn)
    rdata = readr[0]
    readha = read2d(haalln)
    hadata = readha[0]
    hahead = readha[1]
    hadatanew = np.subtract(hadata, (np.multiply(rdata, rharatio)))
    fits.writeto(haallrsub, hadatanew, header=hahead)

    print("Calculating Ha sky:")
    os.system('dmcopy "' + haallrsub + '[(#1,#2)=region(' + skyvalfilen + ')]" ' + haallrsubmaskedn)
    harsubinvlist = reader(haallrsubmaskedn)
    harsubinvlistmod = removezeros(harsubinvlist)
    haskyval = np.median(harsubinvlistmod)
    haskyerr = np.std(harsubinvlistmod)
    print("DONE!")

    readsii = read2d(siialln)
    siidata = readsii[0]
    siihead = readsii[1]
    siidatanew = np.subtract(siidata, (np.multiply(rdata, rsiiratio)))
    fits.writeto(siiallrsub, siidatanew, header=siihead)

    print("Calculating SII sky:")
    os.system('dmcopy "' + siiallrsub + '[(#1,#2)=region(' + skyvalfilen + ')]" ' + siiallrsubmaskedn)
    siirsubinvlist = reader(siiallrsubmaskedn)
    siirsubinvlistmod = removezeros(siirsubinvlist)
    siiskyval = np.median(siirsubinvlistmod)
    siiskyerr = np.std(siirsubinvlistmod)
    print("DONE!")

    print("HaSky:%.2f HaSkyErr:%.2f SIISky:%.2f SIISkyErr:%.2f RSky:%.2f RSkyErr:%.2f" % (haskyval, haskyerr, siiskyval, siiskyerr, rskyval, rskyerr))

    print("Masking and Contouring Allstarred images:")
    barcounter = 0
    maxbarcounter = (len(allimages))
    allmcont = progressbar.ProgressBar(maxval=maxbarcounter).start()
    for i in allimages:
        result1 = "intermediate/" + i.strip(".fit") + "_masked.fit"
        result2 = "intermediate/" + i.strip(".fit") + "_contoured.fit"
        result = "intermediate/" + i.strip(".fit") + "_contoured_binned6.fit"
        togoresult = "intermediate/togo/" + i.strip(".fit") + "_contoured_binned6.fit"
        os.system('dmcopy "' + i + '[exclude (#1,#2)=region(' + problemfilen + ')]" ' + result1)
        os.system('dmcopy "' + result1 + '[(#1,#2)=region(' + contourfilen + ')]" ' + result2)
        filer = fits.open(result2)
        naxis1 = str(filer[0].header['NAXIS1'])
        naxis2 = str(filer[0].header['NAXIS2'])
        filer.close()
        os.system('dmregrid ' + result2 + ' ' + result + ' "1:' + naxis1 + ':' + str(binsize) + ',1:' + naxis2 + ':' + str(binsize) + '" rotangle=0 rotxcenter=0 rotycenter=0 xoffset=0 yoffset=0 npts=0')
        os.system("cp %s %s" % (result, togoresult))
        barcounter += 1
        allmcont.update(barcounter)
    print('\n')

    print("Masking and Contouring PreAllstar images:")
    barcounter = 0
    maxbarcounter = (len(preimages))
    premcont = progressbar.ProgressBar(maxval=maxbarcounter).start()
    for i in preimages:
        result1 = "intermediate/" + i.strip(".fit") + "_masked.fit"
        result2 = "intermediate/" + i.strip(".fit") + "_contoured.fit"
        result = "intermediate/" + i.strip(".fit") + "_contoured_binned6.fit"
        togoresult = "intermediate/togo/" + i.strip(".fit") + "_contoured_binned6.fit"
        os.system('dmcopy "' + i + '[exclude (#1,#2)=region(' + problemfilen + ')]" ' + result1)
        os.system('dmcopy "' + result1 + '[(#1,#2)=region(' + contourfilen + ')]" ' + result2)
        filer = fits.open(result2)
        naxis1 = str(filer[0].header['NAXIS1'])
        naxis2 = str(filer[0].header['NAXIS2'])
        filer.close()
        os.system('dmregrid ' + result2 + ' ' + result + ' "1:' + naxis1 + ':' + str(binsize) + ',1:' + naxis2 + ':' + str(binsize) + '" rotangle=0 rotxcenter=0 rotycenter=0 xoffset=0 yoffset=0 npts=0')
        os.system("cp %s %s" % (result, togoresult))
        barcounter += 1
        premcont.update(barcounter)
    print('\n')

    print("Masking Allstar images:")
    barcounter = 0
    maxbarcounter = (len(allimages))
    allm = progressbar.ProgressBar(maxval=maxbarcounter).start()
    for i in allimages:
        result1 = "intermediate/" + i.strip(".fit") + "_masked.fit"
        result = "intermediate/" + i.strip(".fit") + "_full_binned6.fit"
        togoresult = "intermediate/togo/" + i.strip(".fit") + "_full_binned6.fit"
        filer = fits.open(i)
        naxis1 = str(filer[0].header['NAXIS1'])
        naxis2 = str(filer[0].header['NAXIS2'])
        filer.close()
        os.system('dmregrid ' + result1 + ' ' + result + ' "1:' + naxis1 + ':' + str(binsize) + ',1:' + naxis2 + ':' + str(binsize) + '" rotangle=0 rotxcenter=0 rotycenter=0 xoffset=0 yoffset=0 npts=0')
        os.system("cp %s %s" % (result, togoresult))
        barcounter += 1
        allm.update(barcounter)
    print('\n')

    print("Masking PreAllstar images:")
    barcounter = 0
    maxbarcounter = (len(preimages))
    prem = progressbar.ProgressBar(maxval=maxbarcounter).start()
    for i in preimages:
        result1 = "intermediate/" + i.strip(".fit") + "_masked.fit"
        result = "intermediate/" + i.strip(".fit") + "_full_binned6.fit"
        togoresult = "intermediate/togo/" + i.strip(".fit") + "_full_binned6.fit"
        filer = fits.open(i)
        naxis1 = str(filer[0].header['NAXIS1'])
        naxis2 = str(filer[0].header['NAXIS2'])
        filer.close()
        os.system('dmregrid ' + result1 + ' ' + result + ' "1:' + naxis1 + ':' + str(binsize) + ',1:' + naxis2 + ':' + str(binsize) + '" rotangle=0 rotxcenter=0 rotycenter=0 xoffset=0 yoffset=0 npts=0')
        os.system("cp %s %s" % (result, togoresult))
        barcounter += 1
        prem.update(barcounter)
    print('\n')

    newdir = imdir + '/intermediate/togo'
    os.chdir(newdir)
    R = [ralln.strip(".fit") + '_contoured_binned6.fit', ralln.strip(".fit") + '_full_binned6.fit', rpren.strip(".fit") + '_contoured_binned6.fit', rpren.strip(".fit") + '_full_binned6.fit']
    Ha = [haalln.strip('.fit') + '_contoured_binned6.fit', haalln.strip('.fit') + '_full_binned6.fit', hapren.strip('.fit') + '_contoured_binned6.fit', hapren.strip('.fit') + '_full_binned6.fit']
    SII = [siialln.strip('.fit') + '_contoured_binned6.fit', siialln.strip('.fit') + '_full_binned6.fit', siipren.strip('.fit') + '_contoured_binned6.fit', siipren.strip('.fit') + '_full_binned6.fit']
    Ha_fins = ['I1a1a.fit', 'I1a2a.fit', 'I1b1a.fit', 'I1b2a.fit']
    Ha_ferrs = ['I1a1b.fit', 'I1a2b.fit', 'I1b1b.fit', 'I1b2b.fit']
    SII_fins = ['I2a1a.fit', 'I2a2a.fit', 'I2b1a.fit', 'I2b2a.fit']
    SII_ferrs = ['I2a1b.fit', 'I2a2b.fit', 'I2b1b.fit', 'I2b2b.fit']

    for i in range(0, len(R)):
        haread = read2d(Ha[i])
        hadata = haread[0]
        hahead = haread[1]
        rread = read2d(R[i])
        rdata = rread[0]
        rhead = rread[1]
        siiread = read2d(SII[i])
        siidata = siiread[0]
        siihead = siiread[1]
        finalhadata = hadata - (rharatio * rdata) + (rharatio * binsize * binsize * rskyval) - (
        binsize * binsize * haskyval)
        np.place(finalhadata, (np.absolute(finalhadata - ((rharatio * binsize * binsize * rskyval) - (binsize * binsize * haskyval)))) < 0.9, 0.0)
        finalhaerr = np.sqrt(hadata + (rharatio * rharatio * rdata) + (rharatio * rharatio * binsize * binsize * binsize * binsize * rskyerr * rskyerr) + (binsize * binsize * binsize * binsize * haskyerr * haskyerr))
        finalsiidata = siidata - (rsiiratio * rdata) + (rsiiratio * binsize * binsize * rskyval) - (
        binsize * binsize * siiskyval)
        np.place(finalsiidata, (np.absolute((rsiiratio * binsize * binsize * rskyval) - (binsize * binsize * siiskyval) - finalsiidata)) < 0.9, 0.0)
        finalsiierr = np.sqrt(siidata + (rsiiratio * rsiiratio * rdata) + (rsiiratio * rsiiratio * binsize * binsize * binsize * binsize * rskyerr * rskyerr) + (binsize * binsize * binsize * binsize * siiskyerr * siiskyerr))
        fits.writeto(Ha_fins[i], finalhadata, header=hahead)
        fits.writeto(SII_fins[i], finalsiidata, header=siihead)
        fits.writeto(Ha_ferrs[i], finalhaerr)
        fits.writeto(SII_ferrs[i], finalsiierr)

    os.chdir(imdir)
    localtime = time.asctime(time.localtime(time.time()))
    localtime = localtime.replace(" ", "_")
    logname = localtime + ".log"
    logfile = open(logname, 'w')
    logfile.write("Ha/R=%.2f SII/R=%.2f binsize=%dx%d imdir=%s" % (
    rharatio, rsiiratio, binsize, binsize, imdir))
    logfile.write("HaSky:%.2f HaSkyErr:%.2f SIISky:%.2f SIISkyErr:%.2f RSky:%.2f RSkyErr:%.2f\n" % (haskyval, haskyerr, siiskyval, siiskyerr, rskyval, rskyerr))
    logfile.close()