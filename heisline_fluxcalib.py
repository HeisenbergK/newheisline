from scipy.integrate import simps
import matplotlib.pyplot as plt
from pylab import annotate
from pandas import *
from scipy.stats import linregress
import numpy as np
from astropy.io import fits
from photutils.psf import IterativelySubtractedPSFPhotometry
from photutils import MMMBackground
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.detection import DAOStarFinder
from astropy.modeling.fitting import LevMarLSQFitter
from zscale import *

heislineversion = 3.0
date = "February 18 2018"

# light speed in angstroems per second
c = 2.998E+18


def setstandardfiles():
    # create a list of standard star names and files
    starcount = input("How many stars?\t")
    counter = 1
    starname = []
    while True:
        if counter > starcount:
            break
        curstar = raw_input("Star number %d:\t" % counter)
        starname.append(curstar)
        counter += 1

    return starname


def getstandardmags(filtname, filtfilename, starname, specfilename):
    # integrate and write file
    f = open("result", 'a')
    f.write("Filter\tStar\tFr. Flux\tW/l Flux\tMagnitude\n")
    f.write("Name\tName\terg/s/cm^2/Hz\terg/s/cm^2/A\tAB\n")
    magmaster = DataFrame(np.zeros([len(filtname), len(starname)]), index=filtname, columns=starname)
    magmastern = DataFrame(np.zeros([len(filtname), len(starname)]), index=filtname, columns=starname)
    for i in range(0, len(filtname)):
        for j in range(0, len(starname)):
            filtold = np.genfromtxt(filtfilename[i], dtype=float, delimiter=',')
            if filtold[0, 0] > filtold[1, 0]:
                filt = np.flipud(filtold)
            else:
                filt = filtold
            star = np.genfromtxt(specfilename[j], dtype=float, delimiter=',')
            plt.clf()
            plt.plot(star[:, 0], star[:, 1], 'g')
            plt.title(str(starname[j]) + " spectrum")
            plt.xlabel(r'Wavelength in $\AA$')
            plt.ylabel(r'Hayes Latham mag (blame massey)')
            plt.show()
            ctable = [-0.163, -0.044, 0.055, 0.309]
            atable = np.subtract(57.5, (np.subtract(8.9, ctable)))
            while True:
                filtband = raw_input("Which band does the %s filter belong in? (B,V,R,I)\t" % str(filtname[i]))
                if filtband == "B" or filtband == "b":
                    filtband = 0
                    break
                elif filtband == "V" or filtband == "v":
                    filtband = 1
                    break
                elif filtband == "R" or filtband == "r":
                    filtband = 2
                    break
                elif filtband == "I" or filtband == "i":
                    filtband = 3
                    break
                else:
                    print("Wrong Input")
            atouse = atable[filtband]
            atouse = float(atouse)
            star[:, 1] = np.power(10, (-(star[:, 1]+atouse)/2.5))  # star spectrum in erg/s/cm^2/Hz
            star[:, 1] = c*(np.divide(star[:, 1], (np.square(star[:, 0]))))
            filt_interp = np.interp(star[:, 0], filt[:, 0], filt[:, 1], left=0.0, right=0.0)
            throughput = filt_interp * star[:, 1]
            plt.clf()
            plt.plot(star[:, 0], throughput, 'r', label=("Throughput from filter " + str(filtname[i])))
            plt.plot(star[:, 0], star[:, 1], 'g', label=(str(starname[j]) + " spectrum"))
            plt.legend()
            plt.title(str(starname[j]) + " spectrum and throughput")
            plt.xlabel(r'Wavelength in $\AA$')
            plt.ylabel(r'Energy in $\frac{erg}{s \cdot{} cm^{2} \cdot{} \AA}$')
            plt.show()
            convolved = np.multiply(filt_interp, star[:, 0])
            i1 = simps(throughput, star[:, 0])
            i2 = simps(filt_interp, star[:, 0])
            i3 = simps(convolved, star[:, 0])
            fluxav = i1 / i2  # wavelength flux in erg/s/cm^2/A
            flux = i1  # flux in erg/s/cm^2
            lpivot = i3 / i2  # effective wavelength of filter
            print("lpivot=%s I2=%s Fl=%s Flav=%s" % (str(lpivot), str(i2), str(flux), str(fluxav)))
            ffluxav = (lpivot*lpivot*fluxav)/c  # frequency flux in erg/s/cm^2/Hz
            mabav = -(2.5 * np.log10(ffluxav)) - 48.59
            magmaster[starname[j]][filtname[i]] = ffluxav
            magmastern[starname[j]][filtname[i]] = mabav

            f.write("%s\t%s\t%.8g\t%.8g\t%.8g\n" % (filtname[i], starname[j], ffluxav, fluxav, mabav))
    f.close()
    return magmastern


def getinstrumentfluxes(filtname, starname):
    multiples = DataFrame(np.zeros([len(filtname), len(starname)]), index=filtname, columns=starname)
    for i in filtname:
        for j in starname:
            curmultiples = input("How many times did you observe %s in filter %s?\t" % (j, i))
            curmultiples = int(curmultiples)
            multiples[j][i] = curmultiples
    magmaster = DataFrame(np.zeros([len(filtname), len(starname)]), index=filtname, columns=starname, dtype=object)
    airmasses = DataFrame(np.zeros([len(filtname), len(starname)]), index=filtname, columns=starname, dtype=object)
    for i in filtname:
        for j in starname:
            curmult = multiples[j][i]
            curmult = int(curmult)
            mags = []
            ams = []
            for k in range(0, curmult):
                obsn = k+1
                obsn = int(obsn)
                imagename = j+'_'+i+'-'+str(obsn)+'.fit'
                image = fits.open(imagename)
                readdata = image[0].data

                # fwhm = input("Please enter the fwhm of the sources in observation %i of star %s in filter %s:\t"
                #              % (obsn, j, i))
                fwhm = 7

                readhead = image[0].header
                if "AIRMASS" in readhead:
                    am = readhead["AIRMASS"]
                    am = float(am)
                else:
                    am = raw_input("Enter the airmass of the image %s:\t" % str(imagename))
                    am = float(am)
                if "EXPTIME" in readhead:
                    etime = readhead["EXPTIME"]
                    etime = float(etime)
                elif "EXPOSURE" in readhead:
                    etime = readhead["EXPOSURE"]
                    etime = float(etime)
                else:
                    etime = raw_input("Enter the exposure time of the image %s:\t" % str(imagename))
                    etime = float(etime)

                daogroup = DAOGroup(crit_separation=(3*fwhm))
                mmm_bkg = MMMBackground()
                iraffind = DAOStarFinder(threshold=2.5*mmm_bkg(readdata), fwhm=fwhm)
                fitter = LevMarLSQFitter()
                gaussian_prf = IntegratedGaussianPRF(sigma=(fwhm/2.355))
                gaussian_prf.sigma.fixed = False
                itr_phot_obj = IterativelySubtractedPSFPhotometry(finder=iraffind, group_maker=daogroup,
                                                                  bkg_estimator=mmm_bkg, psf_model=gaussian_prf,
                                                                  fitter=fitter, fitshape=(3*fwhm, 3*fwhm), niters=2)
                phot_results = itr_phot_obj(readdata)
                newtable = phot_results['id', 'x_fit', 'y_fit', 'flux_fit']
                print(newtable)
                print("Please inspect the image for the standard star")
                plt.clf()
                bkgest = mmm_bkg.calc_background(readdata)
                todisp = np.subtract(readdata, bkgest)
                zmin, zmax = zscale_range(todisp)
                plt.imshow(todisp, cmap='gray', aspect=1, origin='lower', vmin=zmin, vmax=zmax)
                for entry in range(0, len(newtable['id'])):
                    xpos = newtable['x_fit'][entry]
                    ypos = newtable['y_fit'][entry]
                    fluxin = newtable['flux_fit'][entry]
                    num = entry + 1
                    print("%s %s %s %s" % (str(num), str(xpos), str(ypos), str(fluxin)))
                    annotate(str(num), xy=(xpos, ypos), xytext=(xpos+0, ypos+fwhm), fontsize=fwhm, xycoords='data',
                             color='green')
                plt.show()
                num = 1
                while True:
                    num = raw_input("Enter the id of the standard star:\t")
                    try:
                        num = int(num)
                        break
                    except ValueError:
                        print("Try again")
                entry = num - 1
                flux = newtable["flux_fit"][entry]
                mag = (2.5 * np.log10(etime)) - (2.5 * np.log10(flux))
                mags.append(mag)
                ams.append(am)
            magmaster[j][i] = mags
            airmasses[j][i] = ams
    return magmaster, airmasses, multiples


def fluxcalibration(filtname, stdmags, insmags, airmasses, stars):
    fluxcal = DataFrame(np.zeros([len(filtname), 2]), index=filtname, columns=['ZP', 'k'])
    stdmags = stdmags.transpose()
    insmags = insmags.transpose()
    airmasses = airmasses.transpose()
    for i in filtname:
        maststds = stdmags[i]
        mastinst = insmags[i]
        mastams = airmasses[i]
        newstds = []
        newinst = []
        newarms = []
        for j in stars:
            curstd = maststds[j]
            curinsts = mastinst[j]
            curams = mastams[j]
            for k in range(0, len(curinsts)):
                curstdn = curstd
                curinstn = curinsts[k]
                curamn = curams[k]
                newstds.append(curstdn)
                newinst.append(curinstn)
                newarms.append(curamn)
        y = np.subtract(newstds, newinst)
        slope, intercept, r_value, p_value, std_err = linregress(newarms, y)
        zp = intercept
        k = slope
        fluxcal['ZP'][i] = zp
        fluxcal['k'][i] = k
    return fluxcal


def fluxcalibrationpivot(filtname, stdmags, insmags, airmasses, stars, katm):
    fluxcal = DataFrame(np.zeros([len(filtname), 2]), index=filtname, columns=['ZP', 'k'])
    stdmags = stdmags.transpose()
    insmags = insmags.transpose()
    airmasses = airmasses.transpose()
    for i in range(0, len(filtname)):
        currk = katm[i]
        print(currk)
        currk = float(currk)
        i = filtname[i]
        maststds = stdmags[i]
        mastinst = insmags[i]
        mastams = airmasses[i]
        newstds = []
        newinst = []
        newarms = []
        ks = []
        for j in stars:
            curstd = maststds[j]
            curinsts = mastinst[j]
            curams = mastams[j]
            for k in range(0, len(curinsts)):
                curstdn = curstd
                curinstn = curinsts[k]
                curamn = curams[k]
                newstds.append(curstdn)
                newinst.append(curinstn)
                newarms.append(curamn)
                ks.append(currk)
        print(newstds)
        print(newinst)
        print(ks)
        print(newarms)
        y = np.subtract(newstds, newinst)
        u = np.multiply(ks, newarms)
        zps = np.add(y, u)
        zp = np.average(zps)
        k = np.average(ks)
        fluxcal['ZP'][i] = zp
        fluxcal['k'][i] = k
    return fluxcal
