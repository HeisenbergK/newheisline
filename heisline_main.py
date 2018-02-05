import curses
from curses import wrapper
from heisline_calls import *
from heisline_fluxcalib import *
from heisline_skyandbin import *
import fileinput
import csv
import progressbar

heislineversion = 2.7
date = "December 31 2017"

'''
filer = open("heisline_version")
versioncont = []
for line in filer:
    versioncont.append(line.strip("\n"))
filer.close()
date = versioncont[1]
'''

bugger = open("heisline_bugs.log")
buggercont = []
for line in bugger:
    buggercont.append(line.strip("\n"))
bugger.close()

bugs = "None"
for i in buggercont:
    if i.startswith(str(heislineversion)):
        bugs = i

if bugs == "None":
    newbug = "None"
else:
    newbug = bugs.split(';')


def main(stdscr):
    # Clear screen
    stdscr.clear()
    curses.init_pair(1, curses.COLOR_RED, curses.COLOR_WHITE)

    stdscr.addstr(0, 0, ("You are using code with Version: " + str(heislineversion) + " (updated on " + date + ")"),
                  curses.A_BOLD)
    if newbug != "None":
        stdscr.addstr(1, 0, "Bug report for current version:", curses.A_STANDOUT)
        if newbug[1] == "OPPERATIONAL":
            stdscr.addstr(2, 0, "The version is operational", curses.A_BLINK)
        if newbug[1] == "NOTOPERATIONAL":
            stdscr.addstr(2, 0, "The version is not operational", curses.A_BLINK)
        stdscr.addstr(3, 0, ("Bug level:" + newbug[4]), curses.A_BOLD)
        stdscr.addstr(4, 0, "Known bugs:", curses.A_BOLD)
        stdscr.addstr(5, 0, newbug[2], curses.A_BOLD)
        stdscr.addstr(6, 0, "Bug effect:", curses.A_BOLD)
        stdscr.addstr(7, 0, newbug[3], curses.A_BOLD)
        stdscr.addstr(8, 0, "Press ENTER to resume code", curses.A_BOLD)
    else:
        stdscr.addstr(1, 0, "Press ENTER to resume code", curses.A_BOLD)
    stdscr.refresh()
    stdscr.getkey()

wrapper(main)

iraf.noao()
iraf.imred()
iraf.ccdred()
iraf.digiphot()
iraf.daophot()
iraf.images()
iraf.tv()

interactivity = masterinteractivitycheck(0)

while True:
    imagedirectory = raw_input("Please insert the directory of the images:\t")
    if os.path.isdir(imagedirectory):
        break
os.chdir(imagedirectory)

iraf.noao.imred.ccdred.setParam('instrument',"ccddb$kpno/camera.dat")

bias = raw_input("Input Bias List:\t")
biaslist = '@' + bias

wantsbias = ask(0)

if wantsbias == 'y':
    combinebias(biaslist)

if wantsbias == 'y':
    if interactivity == 'y':
        interactivity = masterinteractivitycheck(1)

filtnames = choosefilters()

flats = []
for i in range(0, len(filtnames)):
    flats.append("flatlist_" + filtnames[i])
imagelists = []
for i in range(0, len(filtnames)):
    imagelists.append("images_" + filtnames[i])
print("Please make sure you have the following file lists saved in the directory you are using:")
print(flats)
print(imagelists)
dummy = raw_input()

for i in range(0, len(flats)):
    os.system("cp %s %s" % (flats[i], flats[i]+"_b"))
for i in range(0, len(imagelists)):
    os.system("cp %s %s" % (imagelists[i], imagelists[i]+"_b"))
for i in range(0, len(imagelists)):
    os.system("cp %s %s" % (imagelists[i], imagelists[i]+"_b_f"))
for i in range(0, len(imagelists)):
    os.system("cp %s %s" % (imagelists[i], imagelists[i]+"_b_f_sh"))

for i in range(0, len(flats)):
    filer = fileinput.FileInput(flats[i]+"_b", inplace=True)
    for line in filer:
        print(line.replace(".fit", "_b.fit"))
    os.system("grep -v '^$' %s > tmpremoveline" % (flats[i]+"_b"))
    os.system("rm %s" % (flats[i]+"_b"))
    os.system("mv tmpremoveline %s" % (flats[i]+"_b"))
    os.system("rm tmpremoveline")
    filer.close()
for i in range(0, len(imagelists)):
    filer = fileinput.FileInput(imagelists[i]+"_b", inplace=True)
    for line in filer:
        print(line.replace(".fit", "_b.fit"))
    os.system("grep -v '^$' %s > tmpremoveline" % (imagelists[i]+"_b"))
    os.system("rm %s" % (imagelists[i]+"_b"))
    os.system("mv tmpremoveline %s" % (imagelists[i]+"_b"))
    os.system("rm tmpremoveline")
    filer.close()
for i in range(0, len(imagelists)):
    filer = fileinput.FileInput(imagelists[i]+"_b_f", inplace=True)
    for line in filer:
        print(line.replace(".fit", "_b_f.fit"))
    os.system("grep -v '^$' %s > tmpremoveline" % (imagelists[i]+"_b_f"))
    os.system("rm %s" % (imagelists[i]+"_b_f"))
    os.system("mv tmpremoveline %s" % (imagelists[i]+"_b_f"))
    os.system("rm tmpremoveline")
    filer.close()
for i in range(0, len(imagelists)):
    filer = fileinput.FileInput(imagelists[i]+"_b_f_sh", inplace=True)
    for line in filer:
        print(line.replace(".fit", "_b_f_sh.fit"))
    os.system("grep -v '^$' %s > tmpremoveline" % (imagelists[i]+"_b_f_sh"))
    os.system("rm %s" % (imagelists[i]+"_b_f_sh"))
    os.system("mv tmpremoveline %s" % (imagelists[i]+"_b_f_sh"))
    os.system("rm tmpremoveline")
    filer.close()

if wantsbias == 'y':
    print("Bias-Subtracting Flats")
    maxlength = (len(flats))
    flatbiasbar = progressbar.ProgressBar(maxval=maxlength).start()
    for i in range(0, len(flats)):
        subtractbias(imlist=flats[i], biasimage="BIAS.fit")
        flatbiasbar.update(i)
    print("Bias-Subtracting Images")
    maxlength = (len(imagelists))
    imagebiasbar = progressbar.ProgressBar(maxval=maxlength).start()
    for i in range(0, len(imagelists)):
        subtractbias(imlist=imagelists[i], biasimage="BIAS.fit")
        imagebiasbar.update(i)

if wantsbias == 'y':
    if interactivity == 'y':
        masterinteractivitycheck(2)

wantsflat = ask(1)

if wantsflat == 'y':
    for i in range(0, len(flats)):
        combineflats(flats[i], filtnames[i])

if wantsflat == 'y':
    if interactivity == 'y':
        masterinteractivitycheck(3)

if wantsflat == 'y':
    for i in range(0, len(imagelists)):
        flatfielding(imagelists[i], filtnames[i])

if wantsflat == 'y':
    if interactivity == 'y':
        masterinteractivitycheck(4)

imagelists_b_f = []
imagelists_b_f_sh = []
for item in imagelists:
    imagelists_b_f.append(item + "_b_f")
    imagelists_b_f_sh.append(item + "_b_f_sh")
with open("merged_b_f", "wb") as outfile:
    for f in imagelists_b_f:
        with open(f, "rb") as infile:
            outfile.write(infile.read())
with open("merged_b_f_sh", "wb") as outfile:
    for f in imagelists_b_f_sh:
        with open(f, "rb") as infile:
            outfile.write(infile.read())

wantsalign = ask(2)

if wantsalign == 'y':
    if interactivity == 'y':
        masterinteractivitycheck(5)

with open("merged_b_f", "r") as f:
    referenceim = f.readline()
print("Reference Image:\t" + referenceim)

if wantsalign == 'y':
    aligning(referenceim)

basename = raw_input("Please enter the base name of your images:\t")

wantscombine = ask(3)


# Combination
if wantscombine == 'y':
    for i in range(0, len(imagelists_b_f_sh)):
        combining(imagelist=imagelists_b_f_sh[i], filtname=filtnames[i], basename=basename)
# Combination Done
directory = os.getcwd()


# WCS
wantswcs = ask(4)

if wantswcs == 'y':
    for i in filtnames:
        os.system("mkdir %s" % i)
        os.system("cp %s %s" % (basename+"_"+i+".fit", i))
    roughestimatera = raw_input("Enter a rough estimate for the R.A. in HH:MM:SS :\t")
    roughestimatedec = raw_input("Enter a rough estimate for the Dec. in HH:MM:SS :\t")
    directory = os.getcwd()
    for i in filtnames:
        astrometry(directory=directory, filtname=i, basename=basename, roughestimatera=roughestimatera,
                   roughestimatedec=roughestimatedec)
    os.chdir(directory)
# WCS Done


# ALLSTAR
for i in range(0, len(filtnames)):
    # Does the user want allstar?
    while True:
        hewantsallstar = raw_input("Do you want to allstar the %s image?" % (filtnames[i]))
        if hewantsallstar == "N":
            hewantsallstar = "n"
        if hewantsallstar == "Y":
            hewantsallstar = "y"
        if hewantsallstar == "y":
            break
        if hewantsallstar == "n":
            break
    os.chdir(directory + "/" + filtnames[i])
    photalreadyexists = os.path.isfile("photometry1.mag")
    os.chdir(directory)
    if (hewantsallstar == "n") and not photalreadyexists:
        iraf.unlearn("display")
        iraf.unlearn("imexamine")
        os.chdir(directory + "/" + filtnames[i])
        exptime = input("Enter exposure time of the %s image in seconds" % filtnames[i])
        exptime = int(exptime)
        print("Now the %s image will be displayed. Please use the a keystroke on stars across the image to calculate the FWHM, and the m keystroke on background spots to calculate the backgound standard deviation. When you are done press q" %filtnames[i])
        os.system("ds9 &")
        time.sleep(10)
        iraf.images.tv.display(image=(basename+"_"+filtnames[i]+".fit"), frame=1, fill="Yes")
        iraf.images.tv.imexamine(input=(basename+"_"+filtnames[i]+".fit"), frame=1, keeplog="No")
        userpsf = input("Please enter your estimate of the FWHM (based on the above results)\t")
        userbgsigma = input("Please enter your estimate of the backgound standard deviation (based on the above results)\t")
        iraf.unlearn("display")
        iraf.unlearn("datapars")
        iraf.unlearn("findpars")
        iraf.unlearn("daofind")
        iraf.unlearn("centerpars")
        iraf.unlearn("fitskypars")
        iraf.unlearn("photpars")
        iraf.unlearn("phot")
        iraf.unlearn("daopars")
        iraf.unlearn("psf")
        iraf.unlearn("pstselect")
        iraf.unlearn("allstar")
        iraf.unlearn("pfmerge")
        iraf.noao.digiphot.daophot.datapars.setParam('scale', 1)
        iraf.noao.digiphot.daophot.datapars.setParam('fwhmpsf', userpsf)
        iraf.noao.digiphot.daophot.datapars.setParam('sigma', userbgsigma)
        iraf.noao.digiphot.daophot.datapars.setParam('noise', "poisson")
        iraf.noao.digiphot.daophot.datapars.setParam('emission', "Yes")
        iraf.noao.digiphot.daophot.datapars.setParam('readnoise', 8.1)
        iraf.noao.digiphot.daophot.datapars.setParam('epadu', 2.867)
        iraf.noao.digiphot.daophot.datapars.setParam('itime', exptime)
        iraf.noao.digiphot.daophot.findpars.setParam('threshold', 4)
        wdihtdtn = os.getcwd()
        print(wdihtdtn)
        print(basename+"_"+filtnames[i]+"_W.fit")
        iraf.noao.digiphot.daophot.daofind(image=(wdihtdtn+"/"+basename+"_"+filtnames[i]+"_W.fit"),
                                           output="catalog1.coo", boundary="nearest", constant=0, verify="No")
        iraf.noao.digiphot.daophot.centerpars.setParam('calgorithm',"centroid")
        iraf.noao.digiphot.daophot.centerpars.setParam('cbox',5)
        iraf.noao.digiphot.daophot.fitskypars.setParam('salgorithm',"mode")
        iraf.noao.digiphot.daophot.fitskypars.setParam('annulus',(3*userpsf))
        iraf.noao.digiphot.daophot.fitskypars.setParam('dannulus',8)
        iraf.noao.digiphot.daophot.fitskypars.setParam('skyvalue',0)
        iraf.noao.digiphot.daophot.fitskypars.setParam('smaxiter',20)
        iraf.noao.digiphot.daophot.photpars.setParam('weighting',"constant")
        iraf.noao.digiphot.daophot.photpars.setParam('apertures',((2*userpsf)+3))
        iraf.noao.digiphot.daophot.phot(image=(basename+"_"+filtnames[i]+"_W.fit"), coords="catalog1.coo", output="photometry1.mag", verify="No")
    if hewantsallstar == "y":
        iraf.unlearn("display")
        iraf.unlearn("imexamine")
        os.chdir(directory + "/" + filtnames[i])
        while True:
            crowded = raw_input("Do you think the %s image is crowded? y/n\t" % filtnames[i])
            if crowded == "y":
                crowded = "Yes"
                break
            if crowded == "n":
                crowded = "No"
                break
        exptime = input("Enter exposure time of the %s image in seconds" % filtnames[i])
        exptime = int(exptime)
        print("Now the %s image will be displayed. Please use the a keystroke on stars across the image to calculate the FWHM, and the m keystroke on background spots to calculate the backgound standard deviation. When you are done press q" %filtnames[i])
        os.system("ds9 &")
        time.sleep(10)
        iraf.images.tv.display(image=(basename+"_"+filtnames[i]+".fit"), frame=1, fill="Yes")
        iraf.images.tv.imexamine(input=(basename+"_"+filtnames[i]+".fit"), frame=1, keeplog="No")
        userpsf = input("Please enter your estimate of the FWHM (based on the above results)\t")
        userbgsigma = input("Please enter your estimate of the backgound standard deviation (based on the above results)\t")
        iraf.unlearn("display")
        iraf.unlearn("datapars")
        iraf.unlearn("findpars")
        iraf.unlearn("daofind")
        iraf.unlearn("centerpars")
        iraf.unlearn("fitskypars")
        iraf.unlearn("photpars")
        iraf.unlearn("phot")
        iraf.unlearn("daopars")
        iraf.unlearn("psf")
        iraf.unlearn("pstselect")
        iraf.unlearn("allstar")
        iraf.unlearn("pfmerge")
        iraf.noao.digiphot.daophot.datapars.setParam('scale',1)
        iraf.noao.digiphot.daophot.datapars.setParam('fwhmpsf',userpsf)
        iraf.noao.digiphot.daophot.datapars.setParam('sigma',userbgsigma)
        iraf.noao.digiphot.daophot.datapars.setParam('noise',"poisson")
        iraf.noao.digiphot.daophot.datapars.setParam('emission',"Yes")
        iraf.noao.digiphot.daophot.datapars.setParam('readnoise',8.1)
        iraf.noao.digiphot.daophot.datapars.setParam('epadu',2.867)
        iraf.noao.digiphot.daophot.datapars.setParam('itime',exptime)
        iraf.noao.digiphot.daophot.findpars.setParam('threshold',4)
        wdihtdtn=os.getcwd()
        print(wdihtdtn)
        print(basename+"_"+filtnames[i]+"_W.fit")
        iraf.noao.digiphot.daophot.daofind(image=(wdihtdtn+"/"+basename+"_"+filtnames[i]+"_W.fit"), output="catalog1.coo", boundary="nearest", constant=0, verify="No")
        iraf.noao.digiphot.daophot.centerpars.setParam('calgorithm',"centroid")
        iraf.noao.digiphot.daophot.centerpars.setParam('cbox',5)
        iraf.noao.digiphot.daophot.fitskypars.setParam('salgorithm',"mode")
        iraf.noao.digiphot.daophot.fitskypars.setParam('annulus',(3*userpsf))
        iraf.noao.digiphot.daophot.fitskypars.setParam('dannulus',8)
        iraf.noao.digiphot.daophot.fitskypars.setParam('skyvalue',0)
        iraf.noao.digiphot.daophot.fitskypars.setParam('smaxiter',20)
        iraf.noao.digiphot.daophot.photpars.setParam('weighting',"constant")
        iraf.noao.digiphot.daophot.photpars.setParam('apertures',((2*userpsf)+3))
        iraf.noao.digiphot.daophot.phot(image=(basename+"_"+filtnames[i]+"_W.fit"), coords="catalog1.coo", output="photometry1.mag", verify="No")
        notneeded = raw_input("You are now to select 90 stars that are very well-defined, based on their mesh profiles. The mesh profile of the next star is presented by pressing the n key on the ds9 window. The star is accepted by pressing a on the graph window or declined by pressing d on the graph window. When a message on the command window says that the max number of stars is reached, please go to the ds9 window and press q until the code resumes. If you happen to choose a wrong star, panic not. There will be a confirmation step later, where you will be able to remove the star from the list. Keep in mind that you have to note the ID of the star to do so. When you think you are ready press enter")
        os.system("ds9 &")
        time.sleep(10)
        iraf.images.tv.display(image=(basename+"_"+filtnames[i]+".fit"), frame=1, fill="Yes")
        iraf.noao.digiphot.daophot.daopars.setParam('function', "auto")
        iraf.noao.digiphot.daophot.daopars.setParam('varorder', 2)
        iraf.noao.digiphot.daophot.daopars.setParam('nclean', 0)
        iraf.noao.digiphot.daophot.daopars.setParam('saturated', "No")
        iraf.noao.digiphot.daophot.daopars.setParam('matchrad', 3)
        iraf.noao.digiphot.daophot.daopars.setParam('psfrad', ((2*userpsf)+3))
        iraf.noao.digiphot.daophot.daopars.setParam('fitrad', (2*userpsf))
        iraf.noao.digiphot.daophot.daopars.setParam('recenter', "Yes")
        iraf.noao.digiphot.daophot.daopars.setParam('fitsky', "Yes")
        iraf.noao.digiphot.daophot.daopars.setParam('groupsky', crowded)
        iraf.noao.digiphot.daophot.daopars.setParam('sannulus', (3*userpsf))
        iraf.noao.digiphot.daophot.daopars.setParam('wsannulus', 8)
        iraf.noao.digiphot.daophot.daopars.setParam('maxiter', 100)
        iraf.noao.digiphot.daophot.pstselect(image=(basename+"_"+filtnames[i]+".fit"), photfile="photometry1.mag", pstfile="selection.pst", maxnpsf=90, mkstars="Yes", interactive="Yes", verify="No")
        selectioncleanup = raw_input("If you selected a star you think you shouldnt have selected please now find the file named selection.pst inside the filters folder and delete the line(s) beginning with the stars ID(s)")
        iraf.noao.digiphot.daophot.psf(image=(basename+"_"+filtnames[i]+"_W.fit"), photfile=("photometry1.mag"), pstfile="selection.pst", psfimage="model1.psf", opstfile="selectedbypsf1.pst", groupfile="group1.psg", matchbyid="No", interactive="No", mkstars="No", showplots="No", verify="No")
        iraf.noao.digiphot.daophot.allstar(image=(basename+"_"+filtnames[i]+"_W.fit"), photfile=("photometry1.mag"), psfimage="model1.psf", allstarfile="allstar1.als", rejfile="allstarrejects1.arj", subimage=(basename+"_"+filtnames[i]+"_W_allstarsub1.fit"), verify="No")
        # second mandatory run - 2
        iraf.noao.digiphot.daophot.daofind(image=(basename+"_"+filtnames[i]+"_W_allstarsub1.fit"), output="catalog2.coo", boundary="nearest", constant=0, verify="No")
        iraf.noao.digiphot.daophot.phot(image=(basename+"_"+filtnames[i]+"_W_allstarsub1.fit"), coords="catalog2.coo", output="photometry2_1.mag", verify="No")
        iraf.noao.digiphot.daophot.pfmerge(inphotfiles="allstar1.als, photometry2_1.mag", outphotfile="photometry2.mag")
        iraf.noao.digiphot.daophot.allstar(image=(basename+"_"+filtnames[i]+"_W.fit"), photfile="photometry2.mag", psfimage="model1.psf", allstarfile="allstar2.als", rejfile="allstarrejects2.arj", subimage=(basename+"_"+filtnames[i]+"_W_allstarsub2.fit"), verify="No")
        # third mandatory run - 3
        iraf.noao.digiphot.daophot.daofind(image=(basename+"_"+filtnames[i]+"_W_allstarsub2.fit"), output="catalog3.coo", boundary="nearest", constant=0, verify="No")
        iraf.noao.digiphot.daophot.phot(image=(basename+"_"+filtnames[i]+"_W_allstarsub2.fit"), coords="catalog3.coo", output="photometry3_1.mag", verify="No")
        iraf.noao.digiphot.daophot.pfmerge(inphotfiles="allstar2.als, photometry3_1.mag", outphotfile="photometry3.mag")
        iraf.noao.digiphot.daophot.allstar(image=(basename+"_"+filtnames[i]+"_W.fit"), photfile="photometry3.mag", psfimage="model1.psf", allstarfile="allstar3.als", rejfile="allstarrejects3.arj", subimage=(basename+"_"+filtnames[i]+"_W_allstarsub3.fit"), verify="No")
        os.system("cp %s %s" % ((basename+"_"+filtnames[i]+"_W_allstarsub3.fit"), (basename+"_"+filtnames[i]+"_W_allstarfin.fit")))
        allstarversion = 3
        # while True until interrupt until the user does indeed like the result
        while True:
            os.system("ds9 %s" %(basename+"_"+filtnames[i]+"_W_allstarfin.fit"))
            allstargood = raw_input("Did you like the result on allstar subtracted image?")
            if allstargood == "N":
                allstargood = "n"
            if allstargood == "Y":
                allstargood = "y"
            if allstargood == "y":
                break
            nextversion = allstarversion+1
            lastsubtracted = (basename+"_"+filtnames[i]+"_W_allstarsub"+str(allstarversion)+".fit")
            newcatalog = "catalog"+str(nextversion)+".coo"
            newphotprim = "photometry"+str(nextversion)+"_1.mag"
            lastallstar = "allstar"+str(allstarversion)+".als"
            photmergeimputs = lastallstar+", "+newphotprim
            newphot = "photometry"+str(nextversion)+".mag"
            newallstar = "allstar"+str(nextversion)+".als"
            newrejects = "allstarrejects"+str(nextversion)+".arj"
            newsubtracted = basename+"_"+filtnames[i]+"_W_allstarsub"+str(nextversion)+".fit"
            iraf.daofind(image=lastsubtracted, output=newcatalog, boundary="nearest", constant=0, verify="No")
            iraf.phot(image=lastsubtracted, coords=newcatalog, output=newphotprim, verify="No")
            iraf.pfmerge(inphotfiles=photmergeimputs, outphotfile=newphot)
            iraf.allstar(image=(basename+"_"+filtnames[i]+"_W.fit"), photfile=newphot, psfimage="model1.psf", allstarfile=newallstar, rejfile=newrejects, subimage=newsubtracted, verify="No")
            os.system("rm %s" % (basename+"_"+filtnames[i]+"_W_allstarfin.fit"))
            os.system("cp %s %s" % (newsubtracted, (basename+"_"+filtnames[i]+"_W_allstarfin.fit")))
            allstarversion += 1
        # flipper!!!
        iraf.flpr()
iraf.unlearn("display")
iraf.unlearn("datapars")
iraf.unlearn("findpars")
iraf.unlearn("daofind")
iraf.unlearn("centerpars")
iraf.unlearn("fitskypars")
iraf.unlearn("photpars")
iraf.unlearn("phot")
iraf.unlearn("daopars")
iraf.unlearn("psf")
iraf.unlearn("pstselect")
iraf.unlearn("allstar")
iraf.unlearn("pfmerge")
iraf.unlearn("imexamine")
os.chdir(directory)

# ALLSTAR DONE


# Scaling
wantsscaling = ask(5)

if wantsscaling == 'n':
    os.chdir("scaling")
    packs = []
    with open('scalefactors', 'rb') as f:
        reader = csv.reader(f, delimiter=',')
        packs = list(reader)

if wantsscaling == 'y':
    packs, narrows = scaling(directory=directory, filtnames=filtnames, basename=basename)

wantsscalesub = ask(6)

if wantsscalesub == 'y':
    scalesub(packs, basename)
# Scaling done


os.chdir(directory)


# Flux Calibration
wantsflcal = ask(7)

kabs = []
ZP = []

if wantsflcal == 'not yet supported':
    fluxcaldir = raw_input("Enter folder that contains the standard spectra and filter transmittances")
    os.chdir(fluxcaldir)
    transfiles = []
    for i in narrows:
        transfiles.append(i+'.csv')
    starnames = setstandardfiles()
    specnames = []
    for i in starnames:
        specnames.append(i+'.csv')

    stdmags = getstandardmags(filtname=narrows, filtfilename=transfiles, starname=starnames, specfilename=specnames)
    instmags, ams = getinstrumentfluxes(filtname=narrows, starname=starnames)
    calibrationmaster = fluxcalibration(filtname=narrows, stdmags=stdmags, insmags=instmags, airmasses=ams)

if wantsflcal == 'y' or wantsflcal == 'n':
    for i in narrows:
        dummy = input("Enter the absorption coefficient for filter %s:\t" % i)
        kabs.append(dummy)
        dummy = input("Enter the ZeroPoint for filter %s:\t" % i)
        ZP.append(dummy)
# Flux Calbration Done


# Determining Image names and sky regions
print(filtnames)

while True:
    dummy = 0
    hafilname = raw_input("Out of the filters shown above, which one is the Ha (as is shown above)?\t")
    for i in filtnames:
        if i == hafilname:
            dummy += 1
    if dummy == 1:
        break

while True:
    dummy = 0
    siifilname = raw_input("Out of the filters shown above, which one is the SII (as is shown above)?\t")
    for i in filtnames:
        if i == siifilname:
            dummy += 1
    if dummy == 1:
        break

while True:
    dummy = 0
    rfilname = raw_input("Out of the filters shown above, which one is the sdss-r' (as is shown above)?\t")
    for i in filtnames:
        if i == rfilname:
            dummy += 1
    if dummy == 1:
        break

for i in filtnames:
    os.system("cp "+str(i)+"/"+basename+"_"+str(i)+"_W_allstarfin.fit"+" .")


os.system('mkdir results')
filer = open("versioncontrol", 'rw')
filer.write(str(heislineversion))
filer.close()

harawimname = basename + "_" + hafilname + "_W_allstarfin.fit"
siirawimname = basename + "_" + siifilname + "_W_allstarfin.fit"
rrawimname = basename + "_" + rfilname + "_W_allstarfin.fit"

habinned = basename + "_" + hafilname + "_W_allstarfin_binned.fit"
siibinned = basename + "_" + siifilname + "_W_allstarfin_binned.fit"
rbinned = basename + "_" + rfilname + "_W_allstarfin_binned.fit"
rbinnednext = basename + "_" + rfilname + "_W_allstarfin_binned2.fit"

hasciname = basename + "_" + hafilname + "_W_allstarfin_binned_sci.fit"
siisciname = basename + "_" + siifilname + "_W_allstarfin_binned_sci.fit"
haerrname = basename + "_" + hafilname + "_W_allstarfin_binned_err.fit"
siierrname = basename + "_" + siifilname + "_W_allstarfin_binned_err.fit"

os.system('cp ' + basename + "_" + hafilname + "_W_allstarfin.fit " + "results")
os.system('cp ' + basename + "_" + siifilname + "_W_allstarfin.fit " + "results")
os.system('cp ' + basename + "_" + rfilname + "_W_allstarfin.fit " + "results")
os.chdir('results')
mapdir = os.getcwd()

print("Now the Ha image will be opened. Mark some sky regions on the image and save it as a .reg file in the folder where we are now (%s)." % mapdir)
os.system("ds9 %s" % harawimname)
dummy = raw_input("Enter the name of the file you just saved:\t")
os.system("mv %s hasky.reg" % dummy)

print("Now the SII image will be opened. Mark some sky regions on the image and save it as a .reg file in the folder where we are now (%s)" % mapdir)
os.system("ds9 %s" % siirawimname)
dummy = raw_input("Enter the name of the file you just saved:\t")
os.system("mv %s siisky.reg" % dummy)

print("Now the r' image will be opened. Mark some sky regions on the image and save it as a .reg file in the folder where we are now (%s)" % mapdir)
os.system("ds9 %s" % rrawimname)
dummy = raw_input("Enter the name of the file you just saved:\t")
os.system("mv %s rsky.reg" % dummy)


# Determining done

analysissigma = raw_input("Enter the sigma level to be used for the analysis (I suggest at most 0.5 or 1):\t")
analysissigma = float(analysissigma)

siisigma = raw_input("Enter the sigma level to be used for the photometry of the SII image (I suggest at very most 0.5):\t")
siisigma = float(analysissigma)

binsize = raw_input("Enter the bin size in pixels to be used for the analysis (I suggest at least 20):\t")
binsize = int(binsize)

hak = raw_input("Re-Enter the Ha absorption coefficient to be used for the analysis:\t")
hak = float(hak)

hazp = raw_input("Re-Enter the Ha zero-point to be used for the analysis:\t")
hazp = float(hazp)

siik = raw_input("Re-Enter the SII absorption coefficient to be used for the analysis:\t")
siik = float(siik)

siizp = raw_input("Re-Enter the SII zero-point to be used for the analysis:\t")
siizp = float(siizp)

haf = raw_input("Re-Enter the Ha - r' scaling factor (hopefully computed above) to be used for the analysis:\t")
haf = float(haf)

siif = raw_input("Re-Enter the SII - r' scaling factor (hopefully computed above) to be used for the analysis:\t")
siif = float(siif)

# bin and get mag and error for ha and sii
hamag, hamagerr = narrowtot(harawimname, "hasky.reg", rrawimname, "rsky.reg", hak, hazp, analysissigma, binsize, haf,
                            habinned, rbinned, hasciname, haerrname)
siimag, siimagerr = narrowtot(siirawimname, "siisky.reg", rrawimname, "rsky.reg", siik, siizp, siisigma, binsize, siif,
                              siibinned, rbinnednext, siisciname, siierrname)

filer = open("SNRphotometry.cat", 'rw')
filer.write(basename)
filer.write('Ha AB Magnitude\tHa AB Magnitude Uncertainty')
filer.write('%.8f\t%.8f' % (hamag, hamagerr))
filer.write('SII AB Magnitude\tSII AB Magnitude Uncertainty')
filer.write('%.8f\t%.8f' % (siimag, siimagerr))
filer.close()
