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
for file in filenames:
    print("["+str(i)+"/"+str(numfiles)+"]("+str(int(float(i)/float(numfiles)*100))+"%) "+os.path.basename(os.path.dirname(os.path.abspath(file))) + " " + os.path.basename(file))
    current = openPickle(file)
    samples_a1.append([ math.log(x) for x in current[0] ])
    samples_a2.append([ math.log(x) for x in current[1] ])
    samples_b1.append([ math.log(x) for x in current[2] ])
    samples_b2.append([ math.log(x) for x in current[3] ])

    id = os.path.basename(file).split('-')[0]
    ids.append(id)
    redshift = getCatalogValue(catalog,int(id),r'REDSHIFT')
    #redshift = getCatalogValue(catalog,int(id),r'rN')
    # need to decompose redshift = getCatalogValue(catalog,int(id),r'TARGET_FLAG_TARGET')
    #redshift = getCatalogValue(catalog,int(id),r'LOGLBOL')
    if len(redshift) == 0:
        redshift = -1
    else:
        redshift = redshift[0]
    #redshift = getCatalogValue(catalog,int(id),r'rN')
    redshifts.append(redshift)
    #redshifts.append(classy.getClass(redshift))
    #classy.getClass(catalog[id][14]))#14]
    #ug.append(catalog[id][4])
    #iz.append(catalog[id][7])
    i += 1
#print(redshifts)

#fig = corner.corner(samples_a1, labels=plt_labels,
#                    quantiles=[0.16, 0.5, 0.84],
#                    show_titles=True, title_fmt=".2g", title_kwargs={"fontsize": 15})
#fig.savefig(os.path.join(outdir,"alpha1_medians_together.png"))
#plt.close(fig)
#
#fig = corner.corner(samples_a2, labels=plt_labels,
#                    quantiles=[0.16, 0.5, 0.84],
#                    show_titles=True, title_fmt=".2g", title_kwargs={"fontsize": 15})
#fig.savefig(os.path.join(outdir,"alpha2_medians_together.png"))
#plt.close(fig)
#
#fig = corner.corner(samples_b1, labels=plt_labels,
#                    quantiles=[0.16, 0.5, 0.84],
#                    show_titles=True, title_fmt=".2g", title_kwargs={"fontsize": 15})
#fig.savefig(os.path.join(outdir,"beta1_medians_together.png"))
#plt.close(fig)
#
#fig = corner.corner(samples_b2, labels=plt_labels,
#                    quantiles=[0.16, 0.5, 0.84],
#                    show_titles=True, title_fmt=".2g", title_kwargs={"fontsize": 15})
#fig.savefig(os.path.join(outdir,"beta2_medians_together.png"))
#plt.close(fig)

def isDefault(x):
    print(x)
    x = float(x)
    if x == -1 or x == -9.9:
        print("true")
        return True
    return False

fig = plt.figure()
ax = plt.gca()
xs = [ x[0] for x in samples_a1 ]
ys = [ x[1] for x in samples_a1 ]
cs = redshifts
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

cax = ax.scatter( color_x, color_y, c=color_c, alpha = 0.2, edgecolors='none', cmap=jet )
#cax = ax.scatter( xs, ys, c=cs, alpha=0.2, edgecolors='none', cmap=jet )
#cax = ax.scatter( xs, cs, c=ys, alpha=0.2, edgecolors='none', cmap=jet )
cbar = fig.colorbar(cax)
ax.scatter( black_x, black_y, c=black_c, alpha=0.2, edgecolors='none' )

ax.set_title(r'~2200 AGNs and other Variables, CARMA 2-1, log($\alpha_1$) vs.log($\alpha_2$) with redshift')
ax.set_xlabel(r'log($\alpha_1$)')
ax.set_ylabel(r'log($\alpha_2$)')

#ax.set_yscale('log')
#ax.set_xscale('log')
fig.savefig(os.path.join(outdir,"loglog_alpha1_by_2_redshift_medians_together.png"))
plt.show()
plt.close(fig)

#fig = plt.figure()
#ax = plt.gca()
#xs = [ x[0] for x in samples_a1 ]
#ys = [ x[1] for x in samples_a1 ]
#cs = ug
#jet = plt.get_cmap('jet')#'tab20')

#cax = ax.scatter( xs, ys, c=cs, alpha=0.5, edgecolors='none', cmap=jet )
#cbar = fig.colorbar(cax)

#ax.set_title(r'~2200 AGNs and other Variables, CARMA 2-1, log($\alpha_1$) vs.log($\alpha_2$) with U-G Color Difference')
#ax.set_xlabel(r'log($\alpha_1$)')
#ax.set_ylabel(r'log($\alpha_2$)')

#ax.set_yscale('log')
#ax.set_xscale('log')
#fig.savefig(os.path.join(outdir,"loglog_alpha1_by_2_green_medians_together.png"))
#plt.show()
#plt.close(fig)

#fig = plt.figure()
#ax = plt.gca()
#xs = [ x[0] for x in samples_a1 ]
#ys = [ x[1] for x in samples_a1 ]
#cs = iz
#jet = plt.get_cmap('jet')#'tab20')

#cax = ax.scatter( xs, ys, c=cs, alpha=0.5, edgecolors='none', cmap=jet )
#cbar = fig.colorbar(cax)

#ax.set_title(r'~2200 AGNs and other Variables, CARMA 2-1, log($\alpha_1$) vs.log($\alpha_2$) with I-Z Color Difference')
#ax.set_xlabel(r'log($\alpha_1$)')
#ax.set_ylabel(r'log($\alpha_2$)')
#
##ax.set_yscale('log')
#ax.set_xscale('log')
#fig.savefig(os.path.join(outdir,"loglog_alpha1_by_2_iz_medians_together.png"))
#plt.show()
#plt.close(fig)
