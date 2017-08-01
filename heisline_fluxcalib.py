from scipy.integrate import simps
import matplotlib.pyplot as plt
from pandas import *
from scipy.stats import linregress
import numpy as np


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
    magmaster = DataFrame(np.zeros([len(filtname), len(starname)]), index=filtname, columns=starname)
    for i in range(0, len(filtname)):
        for j in range(0, len(starname)):
            filtold = np.genfromtxt(filtfilename[i], dtype=float, delimiter=',')
            if filtold[0, 1] > filtold[1, 1]:
                filt = np.flipud(filtold)
            else:
                filt = filtold
            star = np.genfromtxt(specfilename[j], dtype=float, delimiter=',')
            star[:, 1] = np.power(10, (-(star[:, 1] + 48.59) / (2.5)))
            filt_interp = np.interp(star[:, 0], filt[:, 0], filt[:, 1])
            throughput = filt_interp * star[:, 1]
            i1 = simps(throughput, star[:, 0])
            i2 = simps(filt[:, 1], filt[:, 0])
            i2 = i2
            fluxav = i1 / i2
            flux = i1
            mab = -2.5 * np.log10(flux) - 48.59
            mabav = -2.5 * np.log10(fluxav) - 48.59
            plt.plot(star[:, 0], throughput, 'r')
            plt.plot(star[:, 0], star[:, 1], 'g')
            plt.show()
            magmaster.set_value(filtname[i], starname[j], fluxav)
            f.write("%s %s %.8g %.8g %.8g %.8g\n" % (filtname[i], starname[j], fluxav, flux, mab, mabav))
    f.close()
    return magmaster


def getinstrumentfluxes(filtname, starname):
    magmaster = DataFrame(np.zeros([len(filtname), len(starname)]), index=filtname, columns=starname)
    airmasses = DataFrame(np.zeros([len(filtname), len(starname)]), index=filtname, columns=starname)
    for i in filtname:
        for j in starname:
            magmaster[i][j] = input("Enter the instrumental magnitude (ZP=25) of star %s in filter %s" % (j, i))
            airmasses[i][j] = input("Enter the airmass of star %s in filter %s" % (j, i))
    return magmaster, airmasses


def fluxcalibration(filtname, stdmags, insmags, airmasses):
    print("Guess you used instrumental ZeroPoint = 25")
    fluxcal = DataFrame(np.zeros([len(filtname), 2]), index=filtname, columns=['ZP', 'k'])
    for i in filtname:
        curstds=stdmags[i]
        curinst=insmags[i]
        curams=airmasses[i]
        y=curinst.subtract(curstds, fill_value=0.0)
        y=y.as_matrix()
        slope, intercept, r_value, p_value, std_err = linregress(curams,y)
        ZP = intercept
        k = slope
        fluxcal[i]['ZP'] = ZP
        fluxcal[i]['k'] = k
    return fluxcal
