import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, sproot, splrep, splev
from mpl_toolkits.mplot3d import Axes3D
import pyregion
from matplotlib import axes, cm

heislineversion = 3.1
date = "February 19 2018"


def filterimage():
    return 0

