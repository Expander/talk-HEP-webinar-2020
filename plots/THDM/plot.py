#!/usr/bin/env python
import sys
import numpy as np
import matplotlib as mpl
mpl.use('PDF')
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.patches import Rectangle
import scipy.interpolate
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.ticker import FuncFormatter
from scipy.interpolate import interp1d
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif', weight='normal')
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

def combine_uncert_tuple_2(Mh2L, Mh3L):
         return abs(Mh2L - Mh3L)

def combine_uncert_tuple_10(Mh2L, Mh3L, MhQmin,  MhQmax, MhASmin, MhASmax, MhMTmin, MhMTmax, MhEFT, MhYt):
         return abs(Mh2L - Mh3L) + \
                  max(abs(Mh2L - MhQmin), abs(Mh2L - MhQmax)) + \
                  abs(MhASmax - MhASmin) + \
                  abs(MhMTmax - MhMTmin) + \
                  abs(Mh2L - MhYt) + \
                  abs(Mh2L - MhEFT)

def combine_uncert_tuple(entries, tuple_DMh):
         switcher = {
                  2 : combine_uncert_tuple_2,
                  10: combine_uncert_tuple_10
         }
         return switcher.get(entries, lambda: "nothing")(*tuple_DMh)

def combine_uncert(DMh):
         try:
                  _ = (e for e in DMh)
                  return combine_uncert_tuple(len(DMh), DMh)
         except TypeError:
                  return abs(DMh)

def zip_uncertainties(data, uncerts):
         if (len(uncerts) == 1):
                  return data[:,uncerts[0]]
         elif (len(uncerts) == 2):
                  return zip(data[:,uncerts[0]],
                             data[:,uncerts[1]])
         elif (len(uncerts) == 10):
                  return zip(data[:,uncerts[0]],
                             data[:,uncerts[1]],
                             data[:,uncerts[2]],
                             data[:,uncerts[3]],
                             data[:,uncerts[4]],
                             data[:,uncerts[5]],
                             data[:,uncerts[6]],
                             data[:,uncerts[7]],
                             data[:,uncerts[8]],
                             data[:,uncerts[9]])
         else:
                  print("Error: cannot handle ", len(uncerts), " uncertainties")
                  return []

def draw_contours(filename, xlabel, ylabel, title, cols, outfile, uncerts, xlog=False, ylog=False):
         try:
                  data = np.genfromtxt(filename)
                  # data = data[~np.isnan(data[:,cols[2]])]
         except:
                  print("Error: could not load numerical data from file")
                  return

         x = data[:,cols[0]]
         y = data[:,cols[1]]
         z = data[:,cols[2]]

         if xlog: x = np.log10(x)
         if ylog: y = np.log10(y)

         # uncertainty
         u = list(map(combine_uncert, zip_uncertainties(data, uncerts)))

         N = 100
         xi = np.linspace(x.min(), x.max(), N)
         yi = np.linspace(y.min(), y.max(), N)
         zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
         ui = scipy.interpolate.griddata((x, y), u, (xi[None,:], yi[:,None]), method='cubic')

         fig = plt.figure(figsize=(4.8,4))
         plt.gcf().subplots_adjust(bottom=0.15, left=0.15)
         plt.tick_params(axis='x',direction='in')
         plt.tick_params(axis='y',direction='in')
         plt.gca().yaxis.set_ticks_position('both')
         plt.gca().xaxis.set_ticks_position('both')
         if xlog:
                  plt.gca().set_xticks(np.log10([10**i for i in range(2,6,1)]))
                  plt.gca().set_xticklabels([r'$10^{{{0}}}$'.format(i) for i in range(2,6,1)])
                  # locmin = mpl.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10)) 
                  # plt.gca().xaxis.set_minor_locator(locmin)
                  # plt.gca().xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
                  # plt.axes().xaxis.set_minor_locator(MultipleLocator(10))
         # define the colormap
         cmap = plt.get_cmap("Blues")
         # extract all colors from the .jet map
         cmaplist = [cmap(i) for i in range(cmap.N)]
         # create the new map
         cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
         # define the bins and normalize
         bounds = np.linspace(0,6,13)
         norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
         # plot uncertainty
         plt.imshow(ui, vmin=0, vmax=6, origin='lower',
                    extent=[x.min(), x.max(), y.min(), y.max()],
                    aspect='auto', cmap=cmap,norm=norm)
         clb = plt.colorbar()
         clb.ax.set_ylabel(r'$\Delta M_h/\text{GeV}$',rotation=-90, labelpad=15)
         clb.ax.tick_params(axis='y',direction='in')

         # Optional: use scatter plot with colored points
         # sc = plt.scatter(x,y, s=10, c=u, vmin=min(u), vmax=max(u), lw=0, cmap=plt.cm.get_cmap('RdYlBu'))
         # plt.colorbar(sc)

         csMh = plt.contour(xi, yi, zi,
                            vmin=100, vmax=130, levels=[110,120,125,130,135],
                            cmap=plt.get_cmap("plasma"))
         plt.clabel(csMh, inline=1, fmt=r'%.0f', cmap=plt.get_cmap("plasma"))
         # csDMh = plt.contour(xi, yi, ui, linestyles=['dashed'],
         #                     levels=[0,1,3,5,6],vmin=0,vmax=6, colors='darkblue',alpha=0.7)
         # plt.clabel(csDMh, inline=1, fmt=r'%.0f')

         plt.xlabel(xlabel)
         plt.ylabel(ylabel)
         plt.title(title)

         plt.savefig(outfile,bbox_inches='tight')
         print("saved plot in {}".format(outfile))
         plt.close(fig)

hgthdm_files_MS_MA = [
         r'SplitTHDMTHDMTower_MS_MA_Xt-2.44949_TB-10_Mu-M12-M3-2000',
         r'SplitTHDMTHDMTower_MS_MA_Xt-2.44949_TB-20_Mu-M12-M3-2000'
]

data_dir = r'./'
plot_dir = r'./'

for f in hgthdm_files_MS_MA:
         filename = data_dir + f + '.dat'
         outfile = plot_dir + f + '.pdf'
         draw_contours(filename, r'$M_S/\text{GeV}$', r'$m_A/\text{GeV}$',
                       r'$\tan\beta = 10$, $X_t = \sqrt{6}M_S$, $\mu = M_i = 2\,\text{TeV}$',
                       [0,5,8], outfile, [8,9,10,11,12,13,14,15,16,17],
                       xlog=True, ylog=False)
