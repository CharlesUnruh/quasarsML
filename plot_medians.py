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



#Set this to True for some debug output
debug = False
if len(sys.argv) < 2:
    print("python2.7 do.py file file ...")
    exit()

fit_models = {}
if debug:
    print("[debug] Pricessing "+str(len(sys.argv)-1)+" filenames to determine models needing plotting...")

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
catalog = openPickle("catalogs/catalog.txt.pkl")
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
    redshifts.append(catalog[id][14])#classy.getClass(catalog[id][14]))#14])
    ug.append(catalog[id][4])
    iz.append(catalog[id][7])
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

fig = plt.figure()
ax = plt.gca()
xs = [ x[0] for x in samples_a1 ]
ys = [ x[1] for x in samples_a1 ]
cs = redshifts
jet = plt.get_cmap('jet')#'tab20')

cax = ax.scatter( xs, ys, c=cs, alpha=0.5, edgecolors='none', cmap=jet )
cbar = fig.colorbar(cax)

ax.set_title(r'~2200 AGNs and other Variables, CARMA 2-1, log($\alpha_1$) vs.log($\alpha_2$) with redshift')
ax.set_xlabel(r'log($\alpha_1$)')
ax.set_ylabel(r'log($\alpha_2$)')

#ax.set_yscale('log')
#ax.set_xscale('log')
fig.savefig(os.path.join(outdir,"loglog_alpha1_by_2_redshift_medians_together.png"))
plt.show()
plt.close(fig)

fig = plt.figure()
ax = plt.gca()
xs = [ x[0] for x in samples_a1 ]
ys = [ x[1] for x in samples_a1 ]
cs = ug
jet = plt.get_cmap('jet')#'tab20')

cax = ax.scatter( xs, ys, c=cs, alpha=0.5, edgecolors='none', cmap=jet )
cbar = fig.colorbar(cax)

ax.set_title(r'~2200 AGNs and other Variables, CARMA 2-1, log($\alpha_1$) vs.log($\alpha_2$) with U-G Color Difference')
ax.set_xlabel(r'log($\alpha_1$)')
ax.set_ylabel(r'log($\alpha_2$)')

#ax.set_yscale('log')
#ax.set_xscale('log')
fig.savefig(os.path.join(outdir,"loglog_alpha1_by_2_green_medians_together.png"))
plt.show()
plt.close(fig)

fig = plt.figure()
ax = plt.gca()
xs = [ x[0] for x in samples_a1 ]
ys = [ x[1] for x in samples_a1 ]
cs = iz
jet = plt.get_cmap('jet')#'tab20')

cax = ax.scatter( xs, ys, c=cs, alpha=0.5, edgecolors='none', cmap=jet )
cbar = fig.colorbar(cax)

ax.set_title(r'~2200 AGNs and other Variables, CARMA 2-1, log($\alpha_1$) vs.log($\alpha_2$) with I-Z Color Difference')
ax.set_xlabel(r'log($\alpha_1$)')
ax.set_ylabel(r'log($\alpha_2$)')

#ax.set_yscale('log')
#ax.set_xscale('log')
fig.savefig(os.path.join(outdir,"loglog_alpha1_by_2_iz_medians_together.png"))
plt.show()
plt.close(fig)
