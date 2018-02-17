import numpy
import pointarray


def zscale_range(image_data, contrast=0.25, num_points=600, num_per_row=120):
    """
    Computes the range of pixel values to use when adjusting the contrast
    of FITs images using the zscale algorithm.  The zscale algorithm
    originates in Iraf.  More information about it can be found in the help
    section for DISPLAY in Iraf.
    Briefly, the zscale algorithm uses an evenly distributed subsample of the
    input image instead of a full histogram.  The subsample is sorted by
    intensity and then fitted with an iterative least squares fit algorithm.
    The endpoints of this fit give the range of pixel values to use when
    adjusting the contrast.
    Input:  image_data  -- the array of data contained in the FITs image
                           (must have 2 dimensions)
            contrast    -- the contrast parameter for the zscale algorithm
            num_points  -- the number of points to use when sampling the
                           image data
            num_per_row -- number of points per row when sampling

    Return: 1.) The minimum pixel value to use when adjusting contrast
            2.) The maximum pixel value to use when adjusting contrast
    """

    # check input shape
    if len(image_data.shape) != 2:
        raise ValueError("input data is not an image")

    # check contrast
    if contrast <= 0.0:
        contrast = 1.0

    # check number of points to use is sane
    if num_points > numpy.size(image_data) or num_points < 0:
        num_points = 0.5 * numpy.size(image_data)

    # determine the number of points in each column
    num_per_col = int(float(num_points) / float(num_per_row) + 0.5)

    # integers that determine how to sample the control points
    xsize, ysize = image_data.shape
    row_skip = float(xsize - 1) / float(num_per_row - 1)
    col_skip = float(ysize - 1) / float(num_per_col - 1)

    # create a regular subsampled grid which includes the corners and edges,
    # indexing from 0 to xsize - 1, ysize - 1
    data = []

    for i in xrange(num_per_row):
        x = int(i * row_skip + 0.5)
        for j in xrange(num_per_col):
            y = int(j * col_skip + 0.5)
            data.append(image_data[x, y])

    # actual number of points selected
    num_pixels = len(data)

    # sort the data by intensity
    data.sort()

    # check for a flat distribution of pixels
    data_min = min(data)
    data_max = max(data)
    center_pixel = (num_pixels + 1) / 2

    if data_min == data_max:
        return data_min, data_max

    # compute the median
    if num_pixels % 2 == 0:
        median = data[center_pixel - 1]
    else:
        median = 0.5 * (data[center_pixel - 1] + data[center_pixel])

    # compute an iterative fit to intensity
    pixel_indeces = map(float, xrange(num_pixels))
    points = pointarray.PointArray(pixel_indeces, data, min_err=1.0e-4)
    fit = points.sigmaIterate()

    num_allowed = 0
    for pt in points.allowedPoints():
        num_allowed += 1

    if num_allowed < int(num_pixels / 2.0):
        return data_min, data_max

    # compute the limits
    z1 = median - (center_pixel - 1) * (fit.slope / contrast)
    z2 = median + (num_pixels - center_pixel) * (fit.slope / contrast)

    if z1 > data_min:
        zmin = z1
    else:
        zmin = data_min

    if z2 < data_max:
        zmax = z2
    else:
        zmax = data_max

    # last ditch sanity check
    if zmin >= zmax:
        zmin = data_min
        zmax = data_max

    return zmin, zmax
