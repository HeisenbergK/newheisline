from astropy.io import fits


# Function that will Read the images
def reader(imname):
    imagehdu=fits.open(imname)
    imagedco=imagehdu[0].data
    toreturn=[]
    for i in range(0,len(imagedco)):
        dummy=imagedco[i]
        for j in range(0,len(dummy)):
            toreturn.append(dummy[j])
    imagehdu.close()
    return toreturn


# Function that read images as 2D
def read2d(imname):
    imhdu=fits.open(imname)
    readdata=imhdu[0].data
    readhead=imhdu[0].header
    imhdu.close()
    toreturn=[readdata,readhead]
    return toreturn


# Function that will remove zeros
def removezeros(listname):
    toreturn=[]
    for i in listname:
        if i>0:
            toreturn.append(i)
    return toreturn