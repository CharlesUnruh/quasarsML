#!/usr/bin/python2

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pickle
import seaborn as sns
import os.path
import corner
import copy
import sys
import pandas
from astropy.io import fits

#Set this to True for some debug output
debug = False
if len(sys.argv) < 2:
    print("python2.7 do.py file file ...")
    exit()

fit_models = {}
if debug:
    print("[debug] Pricessing "+str(len(sys.argv)-1)+" filenames to determine models needing plotting...")

def bitFlagString(bitflags):
    remainders = []
    while bitflags != 0:
        remainders.append(bitflags % 2)
        bitflags = bitflags // 2
    string = ""
    for x in remainders[::-1]:
        string += str(x)
    return string

def openCatalog(filename):
    data = fits.open(filename)
    return data[1].data

def getCatalogValue(catalog,ID,field):
    intermediate = catalog[field]
    idx = np.where( catalog['ID']==ID )[0]
    return intermediate[idx]

def openPickle(filename):
    file = filename
    t_pkl = open(file, 'rb')
    params = pickle.load(t_pkl)
    t_pkl.close()
    return params

def samplesPlot(filename):
    coeff = openPickle(filename)
    ndims = coeff.shape[0]
    nwalkers = coeff.shape[1]
    nsteps = coeff.shape[2]
    samples = np.swapaxes(copy.copy( coeff[:, :, nsteps/2:]).reshape((ndims, -1), order='F'), axis1=0, axis2=1)
    return samples

def doubleMADsfromMedian(y,thresh=4):
    # warning: this function does not check for NAs
    # nor does it address issues when
    # more than 50% of your data have identical values
    m = np.median(y)
    abs_dev = np.abs(y - m)
    left_mad = np.median(abs_dev[y <= m])
    right_mad = np.median(abs_dev[y >= m])
    y_mad = left_mad * np.ones(len(y))
    y_mad[y > m] = right_mad
    modified_z_score = 0.6745 * abs_dev / y_mad
    modified_z_score[y == m] = 0
    return np.where(modified_z_score < thresh)[0]

def cleanLocMin(filename):
    samples = samplesPlot(filename)
    for i in range(3):
        #print(len(samples[:,i]))
        keep_indecies = doubleMADsfromMedian(samples[:,i],thresh=3.0)
        samples = samples[keep_indecies,:]
    return samples


catalog = openCatalog('catalogs/combinedCatalog.fits')

outdir = "median_plots"
pwd = os.getcwd()
outdir = os.path.join(pwd,outdir)
try:
    os.stat(outdir)
except:
    os.mkdir(outdir)

filenames = sys.argv[1:]
labels_timescales = [r"$\tau_1$",r"$\tau_2$",r"$\Sigma_{SF}$",r"$\tau_{MA}$"]
labels_chains = [r"$\alpha_1$", r"$\alpha_2$", r"$\beta_0$",r"$\beta_1$"]
i = 0
numfiles = len(filenames)

plt_labels = labels_timescales
if "Chains" in filenames[0]:
    plt_labels = labels_chains

samples_a1 = []
samples_a2 = []
samples_b1 = []
samples_b2 = []
ids = []
redshifts = []


class classifier():
    def __init__(self):
        self.memo = []
    
    def getClass(self,name):
        if name not in self.memo:
            self.memo.append(name)
        return self.memo.index(name)

classy = classifier()

ug = []
iz = []

fwhm_broad_has = []
fwhm_broad_hbs = []
ew_broad_has = []
ew_broad_hbs = []
ew_oiii_5007s = []
ew_fe_mgiis = []
fwhm_civs = []
ew_civs = []
voff_civ_peaks = []

def catalogToValue(x):
    if len(x) == 0:
        return -1
    return x[0]

for file in filenames:
    print("["+str(i)+"/"+str(numfiles)+"]("+str(int(float(i)/float(numfiles)*100))+"%) "+os.path.basename(os.path.dirname(os.path.abspath(file))) + " " + os.path.basename(file))
    current = openPickle(file)
    samples_a1.append([ math.log(x) for x in current[0] ])
    samples_a2.append([ math.log(x) for x in current[1] ])
    samples_b1.append([ math.log(x) for x in current[2] ])
    samples_b2.append([ math.log(x) for x in current[3] ])

    id = os.path.basename(file).split('-')[0]
    ids.append(id)

    fwhm_broad_has.append(catalogToValue(getCatalogValue(catalog,int(id),r'FWHM_BROAD_HA')))
    fwhm_broad_hbs.append(catalogToValue(getCatalogValue(catalog,int(id),r'FWHM_BROAD_HB')))
    ew_broad_has.append(catalogToValue(getCatalogValue(catalog,int(id),r'EW_BROAD_HA')))
    ew_broad_hbs.append(catalogToValue(getCatalogValue(catalog,int(id),r'EW_BROAD_HB')))
    ew_oiii_5007s.append(catalogToValue(getCatalogValue(catalog,int(id),r'EW_OIII_5007')))
    ew_fe_mgiis.append(catalogToValue(getCatalogValue(catalog,int(id),r'EW_FE_MGII')))
    fwhm_civs.append(catalogToValue(getCatalogValue(catalog,int(id),r'FWHM_CIV')))
    ew_civs.append(catalogToValue(getCatalogValue(catalog,int(id),r'EW_CIV')))
    voff_civ_peaks.append(catalogToValue(getCatalogValue(catalog,int(id),r'VOFF_CIV_PEAK')))

    i += 1

def isDefault(x):
    x = float(x)
    if x == -1 or x == -9.9 or x == 0.0:
        return True
    return False

def plot(filename, colors, title):
    fig = plt.figure()
    ax = plt.gca()
    xs = [ x[0] for x in samples_a1 ]
    ys = [ x[1] for x in samples_a1 ]
    cs = colors
    jet = plt.get_cmap('jet')#'tab20')
    #jet = plt.get_cmap('rainbow')#'tab20')
    
    black_x = []
    black_y = []
    black_c = []
    color_x = []
    color_y = []
    color_c = []
    for i in range(len(xs)):
        if isDefault(cs[i]):
            black_x.append(xs[i])
            black_y.append(ys[i])
            black_c.append('black')
        else:
            color_x.append(xs[i])
            color_y.append(ys[i])
            color_c.append(cs[i])
    
    cax = ax.scatter( color_x, color_c, c=color_y, alpha = 0.5, edgecolors='none', cmap=jet )
    #cax = ax.scatter( xs, ys, c=cs, alpha=0.2, edgecolors='none', cmap=jet )
    #cax = ax.scatter( xs, cs, c=ys, alpha=0.2, edgecolors='none', cmap=jet )
    cbar = fig.colorbar(cax)
    #ax.scatter( black_x, black_y, c=black_c, alpha=0.1, edgecolors='none' )
    
    ax.set_title(title)
    ax.set_xlabel(r'log($\alpha_1$)')
    ax.set_ylabel(r'log($\alpha_2$)')
    #ax.set_ylim([-15,25])
    ax.set_xlim([-10,15])
    
    fig.savefig(os.path.join(outdir,filename))#"loglog_alpha1_by_2_redshift_medians_together.png"))
    plt.show()
    plt.close(fig)

title = r'~2200 AGNs and other Variables, CARMA 2-1, log($\alpha_1$) vs.log($\alpha_2$)'
file_pre = "loglog_alpha1_by_2_medians_together_"
plot(file_pre + "fwhm_broad_has.png", fwhm_broad_has, title + r' FWHM Broad HA')
plot(file_pre + "fwhm_broad_bas.png", fwhm_broad_hbs, title + r' FWHM Broad HB')
plot(file_pre + "ew_broad_has.png", ew_broad_has, title + r' EW Broad HA')
plot(file_pre + "ew_broad_hbs.png", ew_broad_hbs, title + r' EW Broad HB')
plot(file_pre + "ew_oiii_5007.png", ew_oiii_5007s, title + r' EW OIII 5007')
plot(file_pre + "ew_fe_mgii.png", ew_fe_mgiis, title + r' EW FE MGII')
plot(file_pre + "fwhm_civ.png", fwhm_civs, title + r' FWHM CIV')
plot(file_pre + "ew_civ.png", ew_civs, title + r' EW CIV')
plot(file_pre + "voff_civ_peak.png", voff_civ_peaks, title + r' VOFF CIV Peak')

