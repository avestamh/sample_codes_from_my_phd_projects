#----------------------------------------------------------------------
# Written by: Sadra Avestan
# Date: Oct 3 2019
# Description: the plot will show the heatmap of the qnn results for the protein secondary structure.
#-----------------------------------------------------------------------
import matplotlib as m
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rc, rcParams
from matplotlib import cm
import matplotlib.patches as patches
import os
#+++++++++++++++++++++++++++++++++
#plt.figure(figsize=(8, 6), dpi=100)
fontsize = 28
fontweight = 'bold'
rc('font', **{'family': 'serif', 'size': fontsize, 'weight' : fontweight, 'serif': ['tgheros']})
rc('text', usetex=True)
fontproperties = {'family':'serif', 'weight' : fontweight, 'size' : fontsize}
#+++++++++++++++++++++++++++++++++
def readdatfile (datfile,geo):
 qntmd = ["cycle","mode","moveindex","move","frame"]
 if geo==0: qnresid = range(1,297)
 if geo==1: qnresid = range(1,319)
 qnlist = qntmd + qnresid
 dfn = pd.read_csv(datfile, delim_whitespace=True, header=None, names=qnlist, skiprows=1, low_memory=False)
# print dfn
 print datfile
 dfn = dfn.drop(['cycle','mode','moveindex','move','frame'], axis=1)
 dfn=dfn.fillna(0)
 if geo==1: 
    print np.arange(1,23)
    dfn=dfn.drop([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22], axis=1)
#   print dfn
    print len(dfn.columns.values)
    dfn.columns = range(1,297)
    print dfn
 dfnss = dfn
 dfnss['s1']=dfnss[range(12,17)].sum(axis=1)
 dfnss['s2']=dfnss[range(19,27)].sum(axis=1)
 dfnss['s3']=dfnss[range(34,38)].sum(axis=1)
 dfnss['h1']=dfnss[range(44,54)].sum(axis=1)
 dfnss['s4']=dfnss[range(60,64)].sum(axis=1)
 dfnss['h2']=dfnss[range(80,95)].sum(axis=1)
 dfnss['s5']=dfnss[range(99,105)].sum(axis=1)
 dfnss['h3']=dfnss[range(105,118)].sum(axis=1)
 dfnss['s6']=dfnss[range(122,129)].sum(axis=1)
 dfnss['h4']=dfnss[range(137,153)].sum(axis=1)
 dfnss['h5']=dfnss[range(156,165)].sum(axis=1)
 dfnss['h6']=dfnss[range(166,176)].sum(axis=1)
 dfnss['h7']=dfnss[range(182,191)].sum(axis=1)
 dfnss['h8']=dfnss[range(192,208)].sum(axis=1)
 dfnss['h9']=dfnss[range(215,231)].sum(axis=1)
 dfnss['s7']=dfnss[range(235,243)].sum(axis=1)
 dfnss['h10']=dfnss[range(248,258)].sum(axis=1)
 dfnss['s8']=dfnss[range(261,270)].sum(axis=1)
 dfnss['h11']=dfnss[range(278,291)].sum(axis=1)
 #print dfnss
 dfnss = dfnss.drop(range(1,297), axis=1)
 #print 'dfn :', dfn
 return dfn,dfnss 
def plotdf (datfile1,datfile2,geo):
 dfn, dfnss=readdatfile(datfile1,geo)
 dfnn, dfnnss=readdatfile(datfile2,geo)
 #print 'dfnss 2nd fun', dfnss
 #print 'transfo :', dfn.transpose()
 #print 'len dfn', len(dfn), 'lentrans', len(dfn.transpose()), 'max dfn',dfn.values.max()
 xmin = 0
 xmax = len(dfnss) 
 ymin = 1
 ymax = len(dfnss.transpose())+1
 print 'ymax', len(dfnss.transpose())
 dfn=dfn.fillna(0)
 dfnn=dfnn.fillna(0)
 print dfn.sum(axis=1)
 print len(dfn), len(dfn.transpose()), dfn.values.max(), dfnn.values.max()
 levels = np.arange(min(dfnss.values.min(),dfnnss.values.min()), max(dfnss.values.max(), dfnnss.values.max()), 10)#, ]dtype=None)
 print 'max', max(dfnss.values.max(), dfnnss.values.max())
 fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True,figsize=(4,6), dpi = 100)
 bottom,top,left,right = 0.11, .92, 0.16, 0.80
 fig.subplots_adjust(bottom=bottom, top=top,left=left, right=right)
 cdict = {
 'red'  :  ( (0.0, 1, 0.0), (0.5, 0.5, 1.), (1., 1, 0.)),
 'green':  ( (0.0, 1., 0.0), (0.2, 0.1, .45), (1., 0.8, .97)),
 'blue' :  ( (0.0, 0.0, 1.0), (0.02, 0.5, .75), (1., 1.0, 1))
}
#  cdict = {'red':   [(0.0,  0.0, 0.0),
#                    (0.5,  1.0, 1.0),
#                    (1.0,  1.0, 1.0)],
#          'green': [(0.0,  0.0, 0.0),
#                    (0.25, 0.0, 0.0),
#                    (0.75, 1.0, 1.0),
#                    (1.0,  1.0, 1.0)],
#          'blue':  [(0.0,  0.0, 0.0),
#                    (0.5,  0.0, 0.0),
#                    (1.0,  1.0, 1.0)]}
 cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
 cmap="gist_rainbow"
# cmap=cm
# print dfnss.columns.values
# print dfnss.index.values
#masked_array = np.ma.array (dfnss, mask=np.isnan(a))
# mask some 'bad' data, in your case you would have: data == 0
 #srfc1=ax1.contourf(dfnss.transpose(),levels,alpha=1,antialiased=True, cmap=cmap, extent=[xmin, xmax, ymin, ymax])#,clim=(cmin, cmax))
 srfc1=ax1.contourf(dfnss.transpose(),levels,alpha=0.9,antialiased=False, cmap=cmap, extent=[xmin, xmax, ymin, ymax])#,clim=(cmin, cmax))
 #scatter_matrix(dfnss, alpha=0.2, figsize=(6, 6), diagonal='kde')
 #ax1.contour(dfnss.transpose(),levels,alpha=1, cmap=cmap,extent=[xmin, xmax, ymin, ymax])#,clim=(cmin, cmax))
 #ax1.contour(dfnss.transpose(),levels,alpha=1, cmap=cmap,extent=[xmin, xmax, ymin, ymax])#,clim=(cmin, cmax))
 #srfc2=ax2.contourf(dfnnss.transpose(),levels,alpha=1,antialiased=True, cmap=cmap,extent=[xmin, xmax, ymin, ymax])#,clim=(cmin, cmax))
 srfc2=ax2.contourf(dfnnss.transpose(),levels,alpha=1,antialiased=False, cmap=cmap,extent=[xmin, xmax, ymin, ymax], extend='both')#,clim=(cmin, cmax))
 #ax2.contour(dfnnss.transpose(),levels,alpha=1, cmap=cmap,extent=[xmin, xmax, ymin, ymax])#,clim=(cmin, cmax))
# ax1.clabel(srfc1, inline=1, fontsize=10)
# print dfnss.transpose()
#fig.subplots_adjust(bottom=0.1, right=0.8, top=0.9)i
#plt.draw()
 p0 = ax1.get_position().get_points().flatten()
 print 'p0', p0
 p1 = ax2.get_position().get_points().flatten()
 print 'p1', p1
#Add an axes at position rect [left, bottom, width, height] where all quantities are in fractions of figure width and height.
#ax_cbar = fig.add_axes([p0[0], 0.97, p1[2]-p0[0], 0.02])
 ax_cbar = fig.add_axes([p1[2]+ 0.02, p1[1] , 0.03, p0[3]-p1[1]])
#cb = fig.colorbar(srfc2, ticks=np.around(levels,decimals=0),format='%1.0f')
 cb =  plt.colorbar(srfc2,cax=ax_cbar,orientation='vertical',ticks=np.around(levels,decimals=0),format='%1.0f') 
#cb=fig.colorbar(srfc1, ax=ax1,ticks=np.around(levels,decimals=0),format='%1.0f'))
#cb=plt.colorbar(srfc1,format='%1.1f')
# aviod transparency in the colorbar
 cb.set_alpha(1)
 cb.draw_all()
# cb.ax.set_title('Number of Contacts', fontsize=14, rotation=270,x=1.0)
 cb.ax.set_ylabel('Number of Contacts', fontsize=14)
# valid keywords are [u'size', u'width', u'color', u'tickdir', u'pad', u'labelsize', u'labelcolor', u'zorder', u'gridOn', u'tick1On', u'tick2On', u'label1On', u'label2On', u'length', u'direction', u'left', u'bottom', u'right', u'top', u'labelleft', u'labelbottom', u'labelright', u'labeltop']
 cb.ax.tick_params(labelsize=10)
# ax2.imshow(dfnn.transpose(), interpolation='sinc', aspect='auto', origin='lower',alpha=.3,cmap=cmap)#
 ax1.set_xlim(0,2500)
 ax2.set_xlim(0,2500)
 ax1.patch.set_facecolor('black')
 ax1.patch.set_alpha(0.5)
 ax2.patch.set_facecolor('black')
 ax2.patch.set_alpha(0.5)
 ax1.set_title(r'Native Contacts', fontsize=18)
 ax2.set_title(r'Non-Native Contacts', fontsize=18)
 ax1.set_ylabel(r'Secondary Structure', fontsize=18, labelpad=10)
 ax2.set_ylabel(r'Secondary Structure', fontsize=18, labelpad=10)
 ax2.set_xlabel(r'Cycles', fontsize=18, labelpad=5)
# ax1.set_xlabel(r'$\theta$ [$^{\circ}$]', fontsize=16, labelpad=10)
 for ax in [ax1,ax2]:
  ax.spines['top'].set_linewidth(1.)
  ax.spines['right'].set_linewidth(1.)
  ax.spines['left'].set_linewidth(1.)
  ax.spines['bottom'].set_linewidth(1.)
# ax.xaxis.set_ticks_position('bottom')
# ax.yaxis.set_ticks_position('left')
# ax.xaxis.set_ticks_position('top')
# majorLocator = MultipleLocator(20)
# ax.yaxis.set_major_locator(majorLocator)
# minorLocator = AutoMinorLocator(4)
# ax.yaxis.set_minor_locator(minorLocator)
  ax.tick_params(axis='y', which='minor', length=2.5,color='black', direction='in', width=0.75)
  ax.tick_params(axis='y', which='major', length=5, color='black', direction='out', width=1.25,pad=0,right='off')
  ax.tick_params(axis='x', which='major', length=5,color='black', direction='in', width=1.25)
  ax.tick_params(axis='x', which='minor', length=2.5, color='black', direction='in', width=0.75,pad=0)
  ax.yaxis.grid(True, 'minor', color='k', linestyle='--',alpha=.2, linewidth=.5)
  ax.xaxis.grid(True, 'major', color='k', linestyle='--',alpha=0.2, linewidth=.5)
# ax.xaxis.grid(True, 'minor', color='white', linestyle='-',alpha=0.5)
# for label in ax.yaxis.get_ticklabels()[::1]:
#    label.set_visible(False)
# for label in ax.yaxis.get_ticklabels()[1::2]:
#    label.set_visible(True)
  majorLocator = MultipleLocator(500)
  ax.xaxis.set_major_locator(majorLocator)
  minorLocator = AutoMinorLocator(5)
  ax.xaxis.set_minor_locator(minorLocator)
# ax.tick_params(which='both', width=1)
# ax.tick_params(which='major', length=6,color='black', direction='out')
# ax.tick_params(which='minor', length=3, color='black', direction='out', pad=0)
  ax.yaxis.set_tick_params(labelsize=10)
  ax.xaxis.set_tick_params(labelsize=14)
# x = np.arange(0,30000/12,5000/12)
# labels = np.arange(0,60,10)
 x = np.arange(0,3000,500)
 labels = np.arange(0,600,100)
 ax2.set_xticks(x)
 ax2.set_xticklabels(labels)
 ax1.set_xticks(x)
#plt.show()
#plt.ylim(1,89)
#fig.tight_layout()
 print np.arange(1,20,1)
#ax2.set_yticks(np.multiply(np.divide([12,16,19,26,34,37,44,53,60,63,80,94,105,117,122,128,137,152,156,164,166,175,182,190,192,207,215,230,235,242,248,257,261,269,278,290],297),20))#, [r"$\mathrm{A^{\prime}}$",'B','C','D','E','F','G'])
#ax2.set_ticks(locs + 0.5, minor=True)
#ax2.set(ticks=locs, ticklabels=labels)
 locs=np.arange(1,20,1)
 labels=[r'$\mathrm{{\ 1}}$',r'$\mathrm{{\ 2}}$',r'$\mathrm{{\ 3}}$',r'$\mathrm{{\ 1}}$',r'$\mathrm{{\ 4}}$',\
 r'$\mathrm{{\ 2}}$',r'$\mathrm{{\ 5}}$',r'$\mathrm{{\ 3}}$',r'$\mathrm{{\ 6}}$',\
 r'$\mathrm{{\ 4}}$',r'$\mathrm{{\ 5}}$',r'$\mathrm{{\ 6}}$',r'$\mathrm{{\ 7}}$',\
 r'$\mathrm{{\ 8}}$',r'$\mathrm{{\ 9}}$',r'$\mathrm{{\ 7}}$',r'$\mathrm{{10}}$',\
 r'$\mathrm{{\ 8}}$',r'$\mathrm{{11}}$']
 ax1.yaxis.set_ticks(locs, minor=True)
 ax1.yaxis.set(ticks=locs+0.5, ticklabels=labels)
 ax2.yaxis.set_ticks(locs, minor=True)
 ax2.yaxis.set(ticks=locs+0.5, ticklabels=labels)
 width=25
 head_width=50
 head_length=0.5
 helixcolor='crimson'
 bbox_props = dict(boxstyle="rarrow,pad=.8", fc="brown", ec="k", lw=0.3)
 bbox_props1 = dict(boxstyle="roundtooth,pad=1.0,tooth_size=.9", fc="mediumblue", ec="k", lw=0.3)
 for ax in [ax1,ax2]:
        x_posi = -260
        ax.text(x_posi, 1.5, "H.", ha="center", va="center", rotation=90, size=2, bbox=bbox_props, alpha=0.1)
        ax.text(x_posi, 2.5, "H.", ha="center", va="center", rotation=90, size=2, bbox=bbox_props, alpha=0.1)
        ax.text(x_posi, 3.5, "H.", ha="center", va="center", rotation=90, size=2, bbox=bbox_props, alpha=0.1)
        ax.text(x_posi, 4.5, "Saaa", ha="center", va="center", rotation=90, size=2, bbox=bbox_props1, alpha=0.1)
        ax.text(x_posi, 5.5, "H.", ha="center", va="center", rotation=90, size=2, bbox=bbox_props, alpha=0.1)
        ax.text(x_posi, 6.5, "Saaa", ha="center", va="center", rotation=90, size=2, bbox=bbox_props1, alpha=0.1)
        ax.text(x_posi, 7.5, "H.", ha="center", va="center", rotation=90, size=2, bbox=bbox_props, alpha=0.1)
        ax.text(x_posi, 8.5, "Saaa", ha="center", va="center", rotation=90, size=2, bbox=bbox_props1, alpha=0.1)
        ax.text(x_posi, 9.5, "H.", ha="center", va="center", rotation=90, size=2, bbox=bbox_props, alpha=0.1)
        ax.text(x_posi, 10.5, "Saa.", ha="center", va="center", rotation=90, size=2, bbox=bbox_props1, alpha=0.1)
        ax.text(x_posi, 11.5, "Saa.", ha="center", va="center", rotation=90, size=2, bbox=bbox_props1, alpha=0.1)
        ax.text(x_posi, 12.5, "Saa.", ha="center", va="center", rotation=90, size=2, bbox=bbox_props1, alpha=0.1)
        ax.text(x_posi, 13.5, "Saa.", ha="center", va="center", rotation=90, size=2, bbox=bbox_props1, alpha=0.1)
        ax.text(x_posi, 14.5, "Saa.", ha="center", va="center", rotation=90, size=2, bbox=bbox_props1, alpha=0.1)
        ax.text(x_posi, 15.5, "Saa.", ha="center", va="center", rotation=90, size=2, bbox=bbox_props1, alpha=0.1)
        ax.text(x_posi, 16.5, "H.", ha="center", va="center", rotation=90, size=2, bbox=bbox_props, alpha=0.1)
        ax.text(x_posi, 17.5, "Saaa", ha="center", va="center", rotation=90, size=2, bbox=bbox_props1, alpha=0.1)
        ax.text(x_posi, 18.5, "H.", ha="center", va="center", rotation=90, size=2, bbox=bbox_props, alpha=0.1)
        ax.text(x_posi, 19.5, "Saaa", ha="center", va="center", rotation=90, size=2, bbox=bbox_props1, alpha=0.1)
 strandc='brown'
 helixcolor='mediumblue'
 colors = [strandc, strandc, strandc,helixcolor,strandc,helixcolor,strandc,helixcolor,strandc,\
 helixcolor,helixcolor,helixcolor,helixcolor,helixcolor,helixcolor,strandc,helixcolor,strandc,helixcolor]
 for xtick, color in zip(ax1.get_yticklabels(), colors):
    xtick.set_color(color)
 for xtick, color in zip(ax2.get_yticklabels(), colors):
    xtick.set_color(color)
# fig.savefig('native-nonnative-time-evolution/%s.png'%traj,dpi=200)
 fig.savefig('/home/yashan/data1/yashan/project3/plots/%s.png'%traj,dpi=200)
analysisdir = "/home/yashan/data1/yashan/project3/analysis/f500/fnter150/temp"
for traj in  sorted(os.listdir(analysisdir)):
# print traj
 if 'geo1' in traj: geo=1; print geo
 if 'geo0' in traj: geo=0; print geo
 qnn = "native-non-native-contacts-coordist-tmd-pull.dat.merged"
 qnnn = "native-non-native-contacts-coordist-nn-tmd-pull.dat.merged"
 datfile1="%s/%s/%s"%(analysisdir,traj,qnn)
 datfile2="%s/%s/%s"%(analysisdir,traj,qnnn)
 qntmd = ["cycle","mode","moveindex","move","frame"]
 if geo==0: qnresid = range(1,297)
 if geo==1: qnresid = range(1,319)
 qnlist = qntmd + qnresid
# print datfile1
# print datfile2
# print "PLOTIGN........." 
 plotdf(datfile1,datfile2,geo)e
