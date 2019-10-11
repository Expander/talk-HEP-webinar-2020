#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import matplotlib.patches as patches

def plot(datafile, outfile, title, label_y, labels, range_y, range_x, curves=[1,1,1], draw_uncertainties = True):
    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

    try:
        dataSS = np.genfromtxt(datafile[0])
        dataFS = np.genfromtxt(datafile[1])
        dataFS2L = np.genfromtxt(datafile[2])
        dataHS = np.genfromtxt('HSSUSY-3L/' + datafile[1])
    except:
        print "Error: could not load numerical data from file"
        exit

    MSSS   = dataSS[:,0]
    MhSS   = dataSS[:,1]
    DMhSS  = dataSS[:,2]
    MSFS   = dataFS[:,0]
    MhFS   = dataFS[:,1]
    DMhFS  = dataFS[:,2]
    MhEFT  = dataHS[:,3]
    DMhEFT = dataHS[:,4]

    MSFS2L  = dataFS2L[:,0]
    MhFS2L  = dataFS2L[:,1]
    DMhFS2L = dataFS2L[:,2]

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', weight='normal')
    fig = plt.figure(figsize=(4,4))
    plt.gcf().subplots_adjust(bottom=0.15, left=0.15) # room for xlabel
    ax = plt.gca()
    ax.set_axisbelow(True)
    ax.xaxis.set_major_formatter(tck.FormatStrFormatter(r'$%g$'))
    ax.yaxis.set_major_formatter(tck.FormatStrFormatter(r'$%g$'))
    ax.get_yaxis().set_tick_params(which='both',direction='in')
    ax.get_xaxis().set_tick_params(which='both',direction='in')
    plt.grid(color='0.5', linestyle=':', linewidth=0.2, dashes=(0.5,1.5))

    plt.xscale('log')
    plt.xlabel(r'$M_S\,/\,\mathrm{GeV}$')
    plt.ylabel(label_y)
    plt.title(title)

    if curves[0]: plt.plot(MSSS, MhSS , 'r-' , linewidth=1.2)
    if curves[1]: plt.plot(MSFS, MhEFT, 'k:' , linewidth=1.2)
    if curves[2]: plt.plot(MSFS2L, MhFS2L, 'b--', linewidth=1.2)
    # plt.plot(MSFS, MhFS , 'b--', linewidth=1.2)

    if draw_uncertainties:
        if curves[0]:
            plt.fill_between(MSSS, MhSS - DMhSS, MhSS + DMhSS,
                             facecolor='red', alpha=0.3, linewidth=0.0)

        if curves[1]:
            plt.fill_between(MSFS, MhEFT - DMhEFT, MhEFT + DMhEFT,
                             facecolor='black', alpha=0.3, linewidth=0.0)

        if curves[2]:
            plt.fill_between(MSFS2L, MhFS2L - DMhFS2L, MhFS2L + DMhFS2L,
                             facecolor='blue', alpha=0.3, linewidth=0.0)

    _labels = []
    if curves[0]: _labels.append(labels[0])
    if curves[1]: _labels.append(labels[1])
    if curves[2]: _labels.append(labels[2])

    leg = plt.legend(_labels,
                     loc='lower right', fontsize=10, fancybox=None, framealpha=None)
    leg.get_frame().set_alpha(1.0)
    leg.get_frame().set_edgecolor('black')
    plt.ylim(range_y)
    plt.xlim(range_x)

    Mhexp = 125.09
    sigma = 0.32
    
    ax.add_patch(
        patches.Rectangle(
            (ax.get_xlim()[0], Mhexp - sigma)  , # (x,y)
            ax.get_xlim()[1] - ax.get_xlim()[0], # width
            2*sigma                            , # height
            color='orange', alpha=0.5, zorder=-1
        )
    )

    plt.text(ax.get_xlim()[0] + 50, Mhexp + sigma + 0.4, r"ATLAS/CMS $\pm1\sigma$", fontsize=8)

    plt.savefig(outfile)
    print "saved plot in ", outfile
    plt.close(fig)

labels = [
    r'\texttt{SOFTSUSY}+\texttt{Himalaya} 3L',
    r'\texttt{FlexibleSUSY}/\texttt{HSSUSY}+\texttt{Himalaya} 3L',
    r'\texttt{FlexibleSUSY}/fixed-order 2L'
]

plot([r'SS_TB-20_Xt--sqrt6.dat', r'FS_TB-20_Xt--sqrt6.dat', r'FS-2L_TB-20_Xt--sqrt6.dat'],
     r'SS_TB-20_Xt--sqrt6.pdf',
     r'$X_t = -\sqrt{6}M_S, \tan\beta = 20$',
     r'$M_h\,/\,\mathrm{GeV}$', labels, [112,135], [500,10000], [1,0,0], True)

plot([r'SS_TB-20_Xt--sqrt6.dat', r'FS_TB-20_Xt--sqrt6.dat', r'FS-2L_TB-20_Xt--sqrt6.dat'],
     r'Mh_MS_TB-20_Xt--sqrt6.pdf',
     r'$X_t = -\sqrt{6}M_S, \tan\beta = 20$',
     r'$M_h\,/\,\mathrm{GeV}$', labels, [112,135], [500,10000], [1,1,0], True)

labels = [
    r'\texttt{FlexibleSUSY}+\texttt{Himalaya} 3L',
    r'',
    r'\texttt{FlexibleSUSY} 2L'
]

plot([r'FS-3L_TB-20_Xt--sqrt6.dat', r'FS_TB-20_Xt--sqrt6.dat', r'FS-2L_TB-20_Xt--sqrt6.dat'],
     r'Mh_2L_vs_3L_MS_TB-20_Xt--sqrt6.pdf',
     r'$X_t = -\sqrt{6}M_S, \tan\beta = 20$',
     r'$M_h\,/\,\mathrm{GeV}$', labels, [112,135], [500,10000], [1,0,1], False)
