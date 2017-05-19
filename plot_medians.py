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

for file in filenames:
    print("["+str(i)+"/"+str(numfiles)+"]("+str(int(float(i)/float(numfiles)*100))+"%) "+os.path.basename(os.path.dirname(os.path.abspath(file))) + " " + os.path.basename(file))
    current = openPickle(file)
    samples_a1.append([ math.log(x) for x in current[0] ])
    samples_a2.append([ math.log(x) for x in current[1] ])
    samples_b1.append([ math.log(x) for x in current[2] ])
    samples_b2.append([ math.log(x) for x in current[3] ])
    i += 1

fig = corner.corner(samples_a1, labels=plt_labels,
                    quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_fmt=".2g", title_kwargs={"fontsize": 15})
fig.savefig(os.path.join(outdir,"alpha1_medians_together.png"))
plt.close(fig)

fig = corner.corner(samples_a2, labels=plt_labels,
                    quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_fmt=".2g", title_kwargs={"fontsize": 15})
fig.savefig(os.path.join(outdir,"alpha2_medians_together.png"))
plt.close(fig)

fig = corner.corner(samples_b1, labels=plt_labels,
                    quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_fmt=".2g", title_kwargs={"fontsize": 15})
fig.savefig(os.path.join(outdir,"beta1_medians_together.png"))
plt.close(fig)

fig = corner.corner(samples_b2, labels=plt_labels,
                    quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_fmt=".2g", title_kwargs={"fontsize": 15})
fig.savefig(os.path.join(outdir,"beta2_medians_together.png"))
plt.close(fig)

fig = plt.figure()
ax = plt.gca()
xs = [ x[0] for x in samples_a1 ]
ys = [ x[1] for x in samples_a1 ]
ax.plot( xs, ys, 'o', c='blue', alpha=0.75, markeredgecolor='none')
#ax.set_yscale('log')
#ax.set_xscale('log')
fig.savefig(os.path.join(outdir,"loglog_alpha1_by_2_medians_together.png"))
plt.show()
plt.close(fig)

##This takes the full paths specified in the arguments and does two things
##First, it parses out the type of model it is, i.e. (1,0) (3,2) etc.
##Second, it stores the full_path as given, in a dict, keyed on the model tuple
#for full_path in sys.argv[1:]:
#    #This assumes UNIX file path delimiters
#    filename = full_path.split('/')[-1]
#    #This also assumes the file-name format will be a certain way:
#    #stuff.otherstuff.modelInfoHere.whatever.could.be.longer
#    segment = filename.split('.')[2]
#    #and within modelInfoHere:
#    #.something-alphas-betasStuff-could-be-longer
#    #Important to note that 'betas' is assumed to take only one digit
#    segment = segment.split('-')
#    alphas = segment[1]
#    #this assumes models' beta value is only one digit long
#    betas = segment[2][0]
#    if debug:
#        print("[debug] "+filename+" ("+alphas+","+betas+")")
#    attrib_tuple = (int(alphas),int(betas))
#    if attrib_tuple not in fit_models:
#        fit_models[attrib_tuple] = []
#    fit_models[attrib_tuple].append(full_path)
#
#if debug:
#    print("[debug] " + "Set of Models = " + str(fit_models.keys))
#
#comparison_pairs = {}
#labels = {}
##This loop constructs two objects used for plotting later
##First is the comparison_pairs map, which is a dict, keyed on the
## model tuple, mapping to a list of 2-lists. The 2-lists are the
## indices as needed to extract the column of data out of the pkl
## files, for the associated attributes, e.g. for a (3,2) model:
## [a0,a1,a2,sigma,b0,b1] gives us indices [0,1,2,4,5] which we pair
## up for later like: [[0,1],[0,2],[0,4],...]
##The second is similar to the first, it creates a dict, keyed on
## the model tuple, that maps to a list of labels, reflecting the
## column headings for the data, used in naming things for the plot.
## for example, (3,2) would produce ["a0","a1","a2","sigma","b0","b1"].
## Constructing it in this way lets us use the same indices for both
## getting our data and labeling it.
#for model in fit_models:
#    alphas = model[0]
#    betas = model[1]
#    #get a list from 0 to alphas
#    alpha_indices = range(alphas)
#    #get a list from alphas+1 to betas+alphas+1
#    #the +1 skips the sigma or SF_inf
#    beta_indices = [ x + alphas+1 for x in range(betas) ]
#    #put the two lists together
#    indices = alpha_indices
#    indices.extend(beta_indices)
#
#    #Do a similar process to the above, but with strings
#    # representing the column headings instead
#    alpha_strings = [ "a"+str(x) for x in range(alphas) ]
#    beta_strings = [ "b"+str(x) for x in range(betas) ]
#    labels[model] = []
#    labels[model].extend(alpha_strings)
#    labels[model].append('sigma')
#    labels[model].extend(beta_strings)
#
#
#    if debug:
#        print("[debug] For model "+str(model))
#        print("[debug] #alphas = " + str(alphas))
#        print("[debug] #betas = " + str(betas))
#        print("[debug] alpha indices = " + str(alpha_indices))
#        print("[debug] beta indices = " + str(beta_indices))
#        print("[debug] indices = " + str(indices))
#        print("[debug] Pairing indices for comparison...")
#
#    #Generate the pairs for comparison, [0,1],[0,2]... etc
#    model_comparison_pairs = []
#    for i in range(len(indices)):
#        for j in range(i+1,len(indices)):
#            model_comparison_pairs.append([indices[i], indices[j]])
#            if debug:
#                print("[debug]("+str(indices[i])+","+str(indices[j])+")")
#    comparison_pairs[model] = model_comparison_pairs
#
##print(comparison_pairs)
#
##The actual plotting, we order by model, then by pair
## This means we have to open each file for every plot,
## but it also means that we don't need to store data
## in memory for long. This might be slower, but will
## work better for a large number of points or quasars.
#for model in fit_models:
#    for pair in comparison_pairs[model]:
#        fig, ax = plt.subplots()
#        colors = iter(cm.rainbow(np.linspace(0, 1, len(fit_models[model]))))
#        plt.title("Model "+str(model)+", "+str(labels[model][pair[0]])+" v "+str(labels[model][pair[1]]))
#        plt.xlabel(str(labels[model][pair[0]]))
#        plt.ylabel(str(labels[model][pair[1]]))
#        axis = [float('Inf'), float('-Inf'), float('Inf'), float('-Inf')]
#        for filename in fit_models[model]:
#            infile = open(filename, 'rb')
#            unpickle = pkl.load(infile)
#            s = [ 1 for x in range(len(unpickle[0])) ]
#            axis[0] = min(axis[0], np.min(unpickle[pair[0]]))
#            axis[1] = max(axis[1], np.max(unpickle[pair[0]]))
#            axis[2] = min(axis[2], np.min(unpickle[pair[1]]))
#            axis[3] = max(axis[3], np.max(unpickle[pair[1]]))
#            plt.scatter(unpickle[pair[0]], unpickle[pair[1]], s=s, color=next(colors), label=filename)
#        #handles, labels = ax.get_legend_handles_labels()
#        #ax.legend(handles, labels)
#
#        #Give the plot 10% of margin, based on the bounding box of all points
#        lrmargin = 0.1 * (axis[1]-axis[0])
#        axis[0] -= lrmargin
#        axis[1] += lrmargin
#        tbmargin = 0.1 * (axis[3]-axis[2])
#        axis[2] -= tbmargin
#        axis[3] += tbmargin
#        plt.axis(axis)
#        plt.savefig("model-"+str(model[0])+"-"+str(model[1])+"-attributes-"+str(labels[model][pair[0]])+"-vs-"+str(labels[model][pair[1]])+".png")
#            
#
#
##
##for arg in sys.argv[3:]:
##    infile = open(arg,'rb')
##    unpickle = pkl.load(infile)
#    
#
#




