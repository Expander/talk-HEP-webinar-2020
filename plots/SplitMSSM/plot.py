#!/usr/bin/env python
import os
import sys
import argparse
import math
import numpy as np
import matplotlib as mpl
mpl.use('PDF')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
#import seaborn
#plt.style.use('seaborn')
#plt.style.use('classic')
#print(plt.style.available)
#plt.style.use('ggplot')
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif', weight='normal')
mpl.rcParams.update({'text.usetex':True, 'text.latex.preamble':[r'\usepackage{amsmath}']})
# mpl.rc('font', size=15)
# For the zoomed inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
#mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
import collections
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import itertools

#print(plt.style.available)
## Paths to the data file
def uncertHSSUSY(data):
    Mh2L = data[:,7]
    Mh3L = data[:,8]
    MhQmin = data[:,9]
    MhQmax = data[:,10]
    MhASmin = data[:,11]
    MhASmax = data[:,12]
    MhMTmin = data[:,13]
    MhMTmax = data[:,14]
    MhEFT = data[:,15]
    MhYt = data[:,16]
    uncert = abs(Mh2L - Mh3L) + np.maximum(abs(Mh2L - MhQmin), abs(Mh2L - MhQmax)) + abs(MhASmax - MhASmin) + abs(MhMTmax - MhMTmin) +  abs(Mh2L - MhYt) + abs(Mh2L - MhEFT)
    return uncert/2.

def uncertSUSYHD(data):
    return data[:,7]
def uncertFH(data):
    return np.array([0.]*len(data[:,6]))

datafiles = {}
datafiles["HSSUSY2"]  = ['SplitMSSM_degenerate_MS_TB-2_Xt-2.44949_Mlow-1500.dat',  4, {"tb":2,"prog":"M3Mi","uncertfunc":uncertHSSUSY}]
# datafiles["HSSUSY10"] = ['SplitMSSM_degenerate_MS_TB-10_Xt-2.44949_Mlow-1500.dat', 4, {"tb":10,"prog":"M3Mi","uncertfunc":uncertHSSUSY}]
# datafiles["HSSUSY20"] = ['SplitMSSM_degenerate_MS_TB-20_Xt-2.44949_Mlow-1500.dat', 4, {"tb":20,"prog":"M3Mi","uncertfunc":uncertHSSUSY}]
datafiles["HSSUSY50"] = ['SplitMSSM_degenerate_MS_TB-50_Xt-2.44949_Mlow-1500.dat', 4, {"tb":50,"prog":"M3Mi","uncertfunc":uncertHSSUSY}]

## Load
def handle_invalid_conv(column):
    if ((column.decode() == 'invalid') or (column.decode() == '-')):
        return -1.
    else:
        return float(column)

def load_data_sets(datafiles):
    datasets = {}
    for iset,(ifile,max_column,iinfo) in datafiles.items():
        handle_invalid = {}
        for i in range(0,max_column):
            handle_invalid[i] = handle_invalid_conv
        print("Loading data from {}".format(ifile))
        # Data as the first item, info dictionary as the second
        datasets[iset] = [np.loadtxt(ifile, converters = handle_invalid), iinfo]
    return datasets
##
def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        fromfile_prefix_chars='@')
    parser.add_argument("--outputfile",
                        default="SplitMSSM_non_degenerate",
                        help="Output plot file name")
    return parser.parse_args()

##
def plot_sets(datasets,styles,outputfile):
    fig, ax = plt.subplots(figsize=(4,4))
    plt.gcf().subplots_adjust(bottom=0.15, left=0.15)
    # scales
    ax.set_xscale('log')
    # labels
    ax.set_ylabel(r'$M_h/\text{GeV}$')
    ax.set_xlabel(r'$M_S/\text{GeV}$')
    ax.set_title(r'split-SUSY, $X_t = \sqrt{6}M_S$')
    # ax limits
    ax.set_ylim(80,160)
    ax.set_xlim(500,10**16)
    ax.set_xticks([10**i for i in np.arange(3,16,2)])
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_tick_params(which='both',direction='in')
    ax.xaxis.set_tick_params(which='both',direction='in')
    # m_h reference line
    mh_line_MS = [500,10**16]
    mh_central = 125.09
    mh_line = np.array([mh_central,mh_central])
    mh_err = math.sqrt(0.21**2+0.11**2)
    ax.plot(mh_line_MS, mh_line, color='orange', linewidth=0)
    ax.fill_between(mh_line_MS,
                    mh_line-mh_err,
                    mh_line+mh_err,
                    facecolor='orange',alpha=1, edgecolor='none')
    # Plot the tb legend
    tb_artist = [
        ax.add_patch(plt.Rectangle((0, 0), 0, 0, facecolor=styles[2], alpha=0.6, linewidth=0.0)),
        ax.add_patch(plt.Rectangle((0, 0), 0, 0, facecolor=styles[50], alpha=0.6, linewidth=0.0))
    ]
    tb_labels = [r"$\tan\beta = 2$",
                 # r"$\tan\beta = 10$",
                 # r"$\tan\beta = 20$",
                 r"$\tan\beta = 50$"]
    leg_tb = ax.legend(tb_artist, tb_labels,
                       loc='lower right', fontsize=8, fancybox=None, framealpha=None)
    plt.gca().add_artist(leg_tb)
    leg_tb.get_frame().set_alpha(1.0)
    leg_tb.get_frame().set_edgecolor('black')

    # Plot the datasets
    for iset in datasets:
        setdata, setinfo = datasets[iset]
        ms = setdata[:,0]
        mh_lower_bound = setdata[:,3]
        mh_upper_bound = setdata[:,4]
        color = styles[setinfo["tb"]]
        linestyle = styles[setinfo["prog"]]
        ax.fill_between(ms[mh_lower_bound > 0],
                        mh_lower_bound[mh_lower_bound > 0],
                        mh_upper_bound[mh_lower_bound > 0],
                        color=color, linewidth=0, alpha=0.5)

    plt.text(10**12, mh_central + 1, r"ATLAS/CMS $\pm1\sigma$", fontsize=6)

    print("Saving {}.pdf".format(outputfile))
    filename = outputfile+".pdf"
    plt.savefig(filename)
    # os.system("pdfcrop {} {}".format(filename,filename))

if __name__ == "__main__":
    cmnd_args = vars(parse_args())
    datasets = load_data_sets(datafiles)
    styles = {
        "M3Mi":{"linestyle":"-"},
        "M32Mi":{"linestyle":"--"},
        "FlexibleEFT":{"linestyle":"dotted"},
        "FH":{"linestyle":"-."},
        2 :"blue",
        10:"deeppink",
        20:"limegreen",
        50:"red"
    }
    plot_sets(datasets,styles,**cmnd_args)
