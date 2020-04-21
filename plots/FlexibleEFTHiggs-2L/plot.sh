#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import matplotlib.patches as patches
import scipy.interpolate

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

def plot(filename, outfile, xlim, ylim, title):
    try:
        data = np.genfromtxt(filename)
    except:
        print "Error: could not load numerical data from file"
        exit

    MS     = data[:,0]
    FEFT1L = data[:,1]
    FEFT2L = data[:,2]
    FO2L   = data[:,3]
    EFT2L  = data[:,4]

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', weight='normal')
    fig = plt.figure(figsize=(4,4))
    plt.gcf().subplots_adjust(bottom=0.15, left=0.15) # room for xlabel
    plt.grid(color='0.5', linestyle=':', linewidth=0.2, dashes=(0.5,1.5))

    ax = plt.gca()
    ax.set_axisbelow(True)
    ax.xaxis.set_major_formatter(tck.FormatStrFormatter(r'$%d$'))
    ax.yaxis.set_major_formatter(tck.FormatStrFormatter(r'$%d$'))
    ax.get_yaxis().set_tick_params(which='both',direction='in')
    ax.get_xaxis().set_tick_params(which='both',direction='in')

    plt.title(title)
    plt.xscale('log')
    plt.xlabel(r'$M_S/\text{GeV}$')
    plt.ylabel(r'$M_h/\text{GeV}$')
    plt.xlim(xlim)
    plt.ylim(ylim)

    # plt.plot(MS, FEFT1L, 'r:' , linewidth=1.2)
    plt.plot(MS, FEFT2L, 'r-' , linewidth=1.2)
    plt.plot(MS, FO2L  , 'b--', linewidth=1.2)
    plt.plot(MS, EFT2L , 'g-.', linewidth=1.2, dashes=(3,2,1,2))

    leg = plt.legend([# r'$\text{FlexibleEFTHiggs/MSSM\ 1L}$',
                      r'$\text{FlexibleEFTHiggs $2\ell$}$',
                      r'$\text{Fixed order $2\ell$}$',
                      r'$\text{EFT $2\ell$}$'],
                     loc='lower right', fontsize=10)
    leg.get_frame().set_alpha(1.0)
    leg.get_frame().set_edgecolor('black')

    Mhexp = 125.10
    sigma = 0.14
    
    ax.add_patch(
        patches.Rectangle(
            (ax.get_xlim()[0], Mhexp - sigma)  , # (x,y)
            ax.get_xlim()[1] - ax.get_xlim()[0], # width
            2*sigma                            , # height
            color='orange', alpha=0.5, zorder=-1
        )
    )

    plt.text(ax.get_xlim()[0] + 30, Mhexp + sigma + 0.4, r"ATLAS/CMS $\pm1\sigma$", fontsize=8)

    plt.tight_layout()
    plt.savefig(outfile)
    print "saved plot in ", outfile
    plt.close(fig)

plot('ms_scan_Xt0_Tb20.dat'     , r'Mh_MS_TB-20_Xt-0.pdf'     ,
     [100, 50000], [80, 131], r'$\tan\beta = 20, X_t = 0$')

plot('ms_scan_Xt-sqrt6_Tb20.dat', r'Mh_MS_TB-20_Xt--sqrt6.pdf',
     [320, 10000], [100, 131], r'$\tan\beta = 20, X_t = -\sqrt{6}M_S$')
