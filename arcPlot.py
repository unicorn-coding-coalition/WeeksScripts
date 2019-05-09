#!/usr/bin/python

# ---------------------------------------------------------------------------------------
# Flexible Arc Plotting code
# Anthony Mustoe
# Weeks lab, UNC
# 2016, 2017, 2018
# 
# Modifications
# Version 2.1 with plot profile option
#
#
#
# ---------------------------------------------------------------------------------------

import sys, os, math, argparse

import RNAtools

import matplotlib 
#matplotlib.use('Agg')
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['xtick.major.width'] = 2.5
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['xtick.minor.size'] = 4 
matplotlib.rcParams['xtick.minor.width'] = 1

import matplotlib.pyplot as plot
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np


class ArcLegend(object):
    """Container for Legend for arc plot"""
    
    def __init__(self, title=None, colors=[], labels=[] ):
        
        assert len(colors)==len(labels), "colors and labels must be the same size" 

        self.title=title
        self.colors = colors
        self.labels = labels
    

    def drawLegend(self, ax, xbound, ybound):
        """Draw legend on axes, computing location based on axes bounds
        ax     = axes object to draw on
        xbound = xboundary of axis
        ybound = yboundary of axis
        """

        spacing = 3
        
        # determine bounding box
        nlines = len(self.colors) + int(self.title is not None) + 1
        height = spacing*nlines       
        width = 15
        if self.title is not None and len(self.title) > width:
            width = len(self.title)
        
        xloc = xbound-width

        if ybound < 0:
            yloc = ybound+height+2
        #    ax.add_patch(patches.Rectangle( (xloc-2, yloc-height+4), width, height, fc='white', lw=1.0))
        else:
            yloc = ybound-height+2
        #    ax.add_patch(patches.Rectangle( (xloc-2, yloc+4), width, height, fc='white', lw=1.0))
        

        # write the title
        if self.title is not None:
            ax.text(xloc, yloc, self.title, horizontalalignment='left',
                    size="10",weight="medium", color="black")
            yloc -= spacing


        for i, c in enumerate(self.colors):
            # self.colors is a list of dict patch specifications, so c is a dict
            ax.add_patch(patches.Rectangle( (xloc, yloc), 1.5, 1.5, **c ) )
            
            # now add label
            ax.text(xloc+2.5, yloc, self.labels[i], horizontalalignment='left',
                    size="8",weight="normal", color="black")
            yloc -= spacing
    
        
        


class ArcPlot(object):

    def __init__(self, title='', fasta=None):
        
        self.title = title

        self.seq = ''
        if fasta:
            self.addFasta(fasta)


        self.seqcolors = ['black']
        self.drawseq = True
        self.grid = False

        self.handleLength = 4*(math.sqrt(2)-1)/3
        self.adjust = 1.4

        self.topPatches = []
        self.botPatches = []
        self.height = [0, 0] # the figure height

        self.intdistance = None

        self.reactprofile = None
        self.reactprofileType = 'SHAPE'

        self.toplegend = None
        self.botlegend = None


    def addFasta(self, fastafile):
        
        read = False
        with open(fastafile) as inp:
            for line in inp:
   
                if line[0]=='>':
                    read = True
                    continue

                line = line.strip()
    
                if read and len(line) > 0:
                    self.seq += line

                                     



    def addArcPath(self, outerPair, innerPair, panel=1, color = 'black', alpha=0.5):
        """ add the arPath object for a given set of parameters"""
        
        outerPair = [outerPair[0]-0.5, outerPair[1]+0.5]
        innerPair = [innerPair[0]+0.5, innerPair[1]-0.5]


        innerRadius = (innerPair[1] - innerPair[0])/2.0
        outerRadius = (outerPair[1] - outerPair[0])/2.0
        
        verts = []
        codes = []

        # outer left
        verts.append( (outerPair[0], 0) )
        codes.append( Path.MOVETO )

        # outer left control 1
        verts.append( (outerPair[0], panel*self.handleLength*outerRadius) )
        codes.append( Path.CURVE4 )

        # outer left control 2
        verts.append( (outerPair[0]+outerRadius*(1-self.handleLength), panel*outerRadius) )
        codes.append( Path.CURVE4 )

        # outer center
        verts.append( (outerPair[0]+outerRadius, panel*outerRadius) )
        codes.append( Path.CURVE4 )

        # outer right control 1
        verts.append( (outerPair[0]+outerRadius*(1+self.handleLength), panel*outerRadius) )
        codes.append( Path.CURVE4 )

        # outer right control 2
        verts.append( (outerPair[1], panel*self.handleLength*outerRadius) )
        codes.append( Path.LINETO )
                       
        # outer right
        verts.append( (outerPair[1], 0) )
        codes.append( Path.LINETO )
                        
        # inner right
        verts.append( (innerPair[1], 0) )
        codes.append( Path.LINETO )
            
        # inner right control 1
        verts.append( (innerPair[1], panel*self.handleLength*innerRadius) )
        codes.append( Path.CURVE4 )

        # inner right control 2
        verts.append( (innerPair[0]+innerRadius*(1+self.handleLength), panel*innerRadius) )
        codes.append( Path.CURVE4 )

        # inner center
        verts.append( (innerPair[0]+innerRadius, panel*innerRadius) )
        codes.append( Path.CURVE4 )

        # inner left control 1 
        verts.append( (innerPair[0]+innerRadius*(1-self.handleLength), panel*innerRadius) )
        codes.append( Path.CURVE4 )

        # inner left control 2
        verts.append( (innerPair[0], panel*self.handleLength*innerRadius) )
        codes.append( Path.LINETO )
 
        # inner left
        verts.append( (innerPair[0], 0) )
        codes.append( Path.LINETO )
 
        # outer left duplicate 
        verts.append( (outerPair[0], 0) )
        codes.append( Path.CLOSEPOLY )
 
        # rescale verts

        if panel == 1:
            adval = self.adjust
        else:
            adval = -(self.adjust-1)
        
        # move arcs away from x-axis
        verts = [ (x,y+adval) for x,y in verts ]
         

        indpath = Path(verts, codes)
        patch = patches.PathPatch(indpath, facecolor=color, alpha=alpha,
                                  linewidth=0, edgecolor='none')
        
        if panel == 1:
            self.topPatches.append(patch)
            if outerRadius > self.height[1]:
                self.height[1] = outerRadius
        
        else:
            self.botPatches.append(patch)
            if outerRadius > self.height[0]:
                self.height[0] = outerRadius



    def plotProfile(self, ax, colthresh = (-10, 0.4, 0.85, 3), heightscale=None):
        """Add a reactivity profile to the axes (self.reactprofile needs to be set)
        
        ax          = axes object to add plot. Expected to be top axis
        colthresh   = thresholds to bin reactivity data. 
                      First element indicates "No-data lower bound"
                      Last element indicates "Max reactivity" -- values are trucated above this
        heightscale = scaling to control height of bars
        """
        
        if self.reactprofile is None:
            return
        
        # check to see if DMS and colthresh defaults not overridden
        if self.reactprofileType == 'DMS' and colthresh == (-10,0.4,0.85,3):
            colthresh = (-10, 0.2, 0.5, 3)
            

        xvals = [ [] for i in range(4) ]
        yvals = [ [] for i in range(4) ]
        
        if heightscale is None:
           heightscale = max(4, min(10, len(self.reactprofile)/50.))

        for x,y in enumerate(self.reactprofile):
            if y is None or y != y or y<colthresh[0]:
                xvals[0].append(x+1)
                yvals[0].append(0.5*heightscale+self.adjust)
            elif y < colthresh[1]:
                xvals[1].append(x+1)
                yvals[1].append(y*heightscale)
            elif y < colthresh[2]:
                xvals[2].append(x+1)
                yvals[2].append(y*heightscale)
            else:
                xvals[3].append(x+1)
                if y > colthresh[3]:
                    yvals[3].append(colthresh[3]*heightscale)
                else:
                    yvals[3].append(y*heightscale)


        ax.bar(xvals[0], yvals[0], alpha=0.7, linewidth=0, color='olive',
               align='center', clip_on=False, bottom=0)
    
        ax.bar(xvals[1], yvals[1], alpha=0.7, linewidth=0, color='black',
               align='center', clip_on=False, bottom=self.adjust)
        ax.bar(xvals[2], yvals[2], alpha=0.7, linewidth=0, color='orange',
               align='center', clip_on=False, bottom=self.adjust)
        ax.bar(xvals[3], yvals[3], alpha=0.7, linewidth=0, color='red',
               align='center', clip_on=False, bottom=self.adjust)
           
        # deal with the axis
        ax.axes.get_yaxis().set_visible(True)
        ax.tick_params(axis='y', direction='out', labelsize=6, left=True, right=False)
        ax.set_yticks( np.array(colthresh[1:])*heightscale+self.adjust )
        
        labels = map(str, colthresh[1:])
        labels[-1] = '>'+labels[-1]
        ax.set_yticklabels( labels )
        ax.set_ylabel('Norm. React.', position=(0,0), ha='left', size=6)
           
        ax.set_frame_on(True)
        for l in ('right','top','bottom'):
            ax.spines[l].set_visible(False)
           
        ax.spines['left'].set_bounds(self.adjust, colthresh[3]*heightscale+self.adjust)


    
    def highlightLowReactivity(self, profile, cutoff=0.0001, window=1):
        """Profile is a np array (0-indexed) of values"""
        
        for i,r in enumerate(profile):
           
            if r<cutoff:
                patch = patches.Rectangle((i+0.5,0),window,4,linewidth=0,fill=True,alpha=0.2)
                self.topPatches.append(patch)     
                patch = patches.Rectangle((i+0.5,-4),window,4,linewidth=0,fill=True,alpha=0.2)
                self.botPatches.append(patch)     
                




    def writePlot(self, outPath="arcs.pdf", bounds=None, write=True, msg=None):
        
        cutbounds = True
        if bounds is None:
            bounds = (0,len(self.seq))
            cutbounds = False
        else:
            bounds = [bounds[0]-0.5, bounds[1]+0.5] # make new copy

        
        doubleplot = len(self.botPatches)>1
        scaleFactor = 0.05
        
        figWidth = (bounds[1]-bounds[0])*scaleFactor
        figHeight = 2*max(self.height)*scaleFactor

        fig = plot.figure( figsize=(figWidth, figHeight)) # 500*scaleFactor
        fig.subplots_adjust(hspace=0.0)

        axT = fig.add_subplot(211)
        axB = None

        for pat in self.topPatches:
            axT.add_patch(pat)
            pat.set_clip_on(cutbounds)

        axT.set_xlim(bounds)
        axT.set_ylim(0, max(self.height))

        axT.set_frame_on(False)
        axT.axes.get_yaxis().set_visible(False)
        
        if self.toplegend:
            self.toplegend.drawLegend(axT, bounds[1], self.height[1])

        if doubleplot:
            axB = fig.add_subplot(212, sharex=axT)
            for pat in self.botPatches:
                axB.add_patch(pat)
                pat.set_clip_on(cutbounds)

            axB.set_xlim(bounds)
            axB.set_ylim(-max(self.height), 0)
            axB.set_frame_on(False)
            axB.axes.get_yaxis().set_visible(False)
            axB.tick_params(axis='both', which='both', top='off',bottom='off', labelbottom='off')
            
            if self.botlegend:
                self.botlegend.drawLegend(axB, bounds[1], -self.height[0]) 



        # format sequence ticks
        axT.get_xaxis().tick_bottom()  
        
        majorLocator = MultipleLocator(100)
        minorLocator = MultipleLocator(20)
        tickinterval = 5

        axT.xaxis.set_major_locator(majorLocator)
        axT.xaxis.set_major_formatter( FormatStrFormatter('%i') )
        axT.xaxis.set_minor_locator(minorLocator)
        axT.xaxis.set_minor_formatter( FormatStrFormatter('%i') )

        majlabels = axT.axes.get_xaxis().get_majorticklabels()
        minlabels = axT.axes.get_xaxis().get_minorticklabels()
        majticks = axT.axes.get_xaxis().get_major_ticks()
        minticks = axT.axes.get_xaxis().get_minor_ticks()
        

        for label in majlabels:
            label.set_size(10)
        
        if bounds[0] == 0: 
            majlabels[1].set_visible(False)
            minlabels[1].set_visible(False)
            majticks[1].set_visible(False)
            minticks[1].set_visible(False)
        
        beginint, r = divmod(bounds[0], 20)
        if r == 0:
            beginint+=1

        for i, label in enumerate(minlabels):
            label.set_size(7)
            if i % tickinterval == beginint:
                label.set_visible(False)


        # plot gridlines
        if self.grid:
            axT.grid(b=True, which="minor", color='black',linestyle='--', alpha=0.5)
            axT.grid(b=True, which="major", color='black',linestyle='-', lw=2, alpha=0.8)

            if doubleplot:
                axB.grid(b=True, which="minor", color='black',linestyle='--', alpha=0.5)
                axB.grid(b=True, which="major", color='black',linestyle='-', lw=2, alpha=0.8)



        # write title (structure names)
        if self.title:
            axT.text(0.0, self.height[1]-2, self.title, horizontalalignment='left',
                      size="24",weight="bold", color="blue")


        # put nuc sequence on axis 
        if self.drawseq:
            fontProp = matplotlib.font_manager.FontProperties(family = "sans-serif", 
                                              style="normal",
                                              weight="extra bold",
                                              size="4")
            
            for i, nuc in enumerate(self.seq):
                
                if i<bounds[0]-1:
                    continue
                elif i>bounds[1]-1:
                    break

                if nuc == "T":
                    nuc = "U"

                try:
                    col = self.seqcolors[i]
                except IndexError:
                    col = "black"
                
                axT.annotate(nuc, xy=(i+0.5, 0),fontproperties=fontProp,color=col,
                             annotation_clip=False,verticalalignment="baseline")
        

        # INDEXING MIGHT BE OFFF -- CHECK!!!
        if self.intdistance is not None:
            xvals = np.arange(bounds[0]+1, bounds[1]+1)
            yvals = self.intdistance[bounds[0]:bounds[1]+1]
            if np.mean(yvals) > 0:
                axT.plot( xvals, yvals, color="black", lw=2)
            else:
                axB.plot( xvals, yvals, color="black", lw=2)

        if self.reactprofile is not None:
            self.plotProfile(axT)

        if msg is not None:
            axT.text(0,1, msg, transform=axT.transAxes)


        # save the figure
        if write:
            fig.savefig(outPath, dpi=100, bbox_inches="tight", transparent=True)
            plot.close(fig)
        else:
            return fig, axT, axB 



    def addCT(self, ctObj, color='black', alpha=0.7, panel=1):
    
        seqlen = len(ctObj.ct)

        if self.seq == '' or self.seq[0] == ' ':
            self.seq = ctObj.seq
        elif len(self.seq) != seqlen:
            print ("Warning:: CT length = {0}; expecting length = {1}").format(seqlen, len(self.seq))
            #self.seq = ctObj.seq
        

        i = 0
        while i < seqlen:
            
            if ctObj.ct[i] > i:

                #outerPair = [i+0.5, ctObj.ct[i]+0.5]
                outerPair = [ i+1, ctObj.ct[i] ]

                # find the right side of helix
                lastPairedNuc = ctObj.ct[i]
                i += 1
                while i<seqlen and (lastPairedNuc-ctObj.ct[i]) == 1:
                    lastPairedNuc = ctObj.ct[i]
                    i += 1
                i -= 1

                #innerPair = [i+1.5, ctObj.ct[i]-0.5] 
                innerPair = [ i+1, ctObj.ct[i] ]

                self.addArcPath(outerPair, innerPair, panel=panel, color=color, alpha=alpha)

            i += 1




    def addDotPlot(self, dpObj, scheme = 'pairprob', panel=1, bins=None):
        """ add dotplot arcs. 
            If scheme is pairprob, plot as pairing probs
            If scheme is corr, plot as correlations
        """
        
        if self.seq == '':
            self.seq = ' '*dpObj.length
        elif len(self.seq) != dpObj.length:
            print ("Warning:: dp file sequence length = {0}; CT/FASTA length = {1}").format(dpObj.length, len(self.seq))
            #raise ValueError('dpObj has different length than expected')


        if scheme=='pairprob':
            refcolors = [ (150,150,150), (255,204,0),  (72,143,205) ,(81, 184, 72) ]
            refalpha = [0.55, 0.65, 0.75, 0.85]

            gradiant = [False, False, False, False]

            if bins is None:
                bins = [0.03, 0.1, 0.3, 0.8, 1.0]
                colors = refcolors
                alpha = refalpha
            else:
                if len(bins) > 5:
                    raise IndexError('{0} PairProb levels specified -- the maximum allowed is 4'.format(len(bins)-1))

                colors = refcolors[-(len(bins)-1):]
                alpha = refalpha[-(len(bins)-1):]


        elif scheme=='corr':
            colors = [(44,123,182), (44,123,182), (171,217,233), (255,255,255), (253,174,97), (215,25,28), (215,25,28)]
            alpha = [1.0, 1.0, 0.0, 0.0, 1.0, 1.0]
            gradiant = [False, True, False, False, True, False]
            
            if bins is None:
                bins = [ -1.0, -0.06, -0.02, 0, 0.02, 0.06, 1.0]
        
        else:
            raise ValueError("Invalid scheme :: {0}".format(scheme))


        # convert colors to percents
        colors = [tuple([y/255.0 for y in x]) for x in colors]
        
        for i in range(len(bins)-1):
            
            # get all pairs within the bin threshold
            if scheme == 'corr':
                cutdp = dpObj.requireProb(10**-bins[i+1], 10**-bins[i])
            else:
                cutdp = dpObj.requireProb(bins[i], bins[i+1])
            cutpairs = cutdp.pairList()
            
            if len(cutpairs) == 0:
                continue
            
            # sort each of the pairs by its strength
            strength = []
            for p in cutpairs:
                elem = cutdp.partfun[p[0]-1]
                pos = ( elem['pair'] == p[1]-1) # make mask
                strength.append(elem['log10'][pos][0]) # get value

            strengthzip = zip(cutpairs, strength)
            strengthzip.sort(key=lambda x: x[1])

            #if scheme == 'pairprob':
            #    strengthzip.reverse()

            cutpairs, strength = zip(*strengthzip)
            
            # iterate through sorted pairs and draw arcs
            for pindex, pair in enumerate(cutpairs):
 
                if gradiant[i]:
                    col = self._colorGrad(strength[pindex], colors[i], colors[i+1], 
                                          bins[i], bins[i+1], log=(scheme=='pairprob'))
                else:
                    col = colors[i]

                #outerPair = (pair[0]-0.5, pair[1]+0.5)
                #innerPair = (pair[0]+0.5, pair[1]-0.5)
                #self.addArcPath(outerPair, innerPair, panel=panel, color=col, alpha=alpha[i])
                self.addArcPath(pair, pair, panel=panel, color=col, alpha=alpha[i])





    def _colorGrad(self, value, colorMin, colorMax, minValue, maxValue, log=False):
        """ return an interpolated rgb color value """
        
        if log:
            value = 10**-value
            minValue = 10**-minValue
            maxValue = 10**-maxValue

        if value > maxValue:
            return colorMax
        elif value < minValue:
            return colorMin
        else:
            v = value - min(maxValue, minValue)
            v /= abs(maxValue-minValue)
        
            col = []
            for i in range(3):
                col.append( v*(colorMax[i] - colorMin[i]) + colorMin[i] )
            
            return col


    def compareCTs(self, refCT, compCT, panel=1):
        
        share = refCT.copy()
        refonly = refCT.copy()
        componly = refCT.copy()
        
        for c in (share, refonly, componly):
            for i in range(len(c.ct)):
                c.ct[i] = 0
    

        for i in xrange(len(refCT.ct)):
            
            if refCT.ct[i] == compCT.ct[i]:
                share.ct[i] = refCT.ct[i]
            else:
                refonly.ct[i] = refCT.ct[i]
                componly.ct[i] = compCT.ct[i]
            

        self.addCT(share, color='green', panel=panel)
        self.addCT(refonly, color='red', panel=panel)
        self.addCT(componly, color='purple', panel=panel)


    def assignSeqColors(self, shapearr):

        if len(shapearr) != len(self.seq):
            raise IndexError("Shape Array does not match length of sequence")

        col = []       
        for x in shapearr:

            if x < -4: 
                col.append( (160,160,160 ) )  # grey
            elif x > 0.85: 
                col.append( (255,0,0) )  # red
            elif 0.85 >= x >0.4: 
                col.append( (255,164,26) ) # orange
            else: 
                col.append( (0,0,0) ) # black

        self.seqcolors = [tuple([y/255.0 for y in x]) for x in col]


    
    def colorSeqByMAP(self, mapfile):
        #hasattr(inp, '__call__')
        
        shape, seq = RNAtools.readSHAPE(mapfile)
        
        if self.seq == '' or self.seq.count(' ')==len(self.seq):
            self.seq = seq
            
        elif len(shape) != len(self.seq):
            raise IndexError("Mapfile does not match length of sequence!")
        
        self.assignSeqColors(shape)
 
    
    def readProfile(self, mapfile):

        shape, seq = RNAtools.readSHAPE(mapfile)
        
        if self.seq == '' or self.seq.count(' ') == len(self.seq):
            self.seq = seq
        
        elif len(shape) != len(self.seq):
            raise IndexError("Mapfile does not match length of sequence!")
        
        self.reactprofile = shape


    def readDMSProfile(self, profilefile):
        
        with open(profilefile) as inp:
            line = inp.readline().split()
            if len(line) == 2:
                ftype=1
            else:
                ftype=2
        

        if ftype==1:
            self.readProfile(profilefile)

        else:
            import plotTools
            profile = plotTools.ReactivityProfile(profilefile)
            profile.normalize(DMS=True)
            self.reactprofile = profile.normprofile
        

        self.reactprofileType = 'DMS'


        
    def old_addInteractionDistance(self, readMatFile, thresh, panel=1):
        
        height = []
        readMat = np.loadtxt(readMatFile)

        for i in xrange(len(self.seq)):
            
            for j in xrange(i,len(self.seq)):
                if readMat[i][j] < thresh:
                    break
            n = float(j-i)

            for j in xrange(i, -1, -1):
                if readMat[i][j] < thresh:
                    break
            m = float(i-j)

            try:
                z = n*math.sin(math.atan(m/n))
            except:
                z = 0

            height.append(panel*(z+self.adjust))

        self.intdistance = height



    def addInteractionDistance(self, readMatFile, thresh, panel=1):
        
        readmat = np.loadtxt(readMatFile)
        
        yvals = np.zeros(len(self.seq))

        for i in xrange(len(self.seq)):

            for j in xrange(i, len(self.seq)):
                if readmat[i,j] >=thresh:

                    idx = int(round((j-i+1)/2))
                    
                    if idx > yvals[i+idx]:
                        yvals[i+idx] = idx
        
        self.intdistance = panel*yvals


#############################################################################




def parseArgs():

    prs = argparse.ArgumentParser()
    prs.add_argument("outputPDF",type=str, help="Name of the output PDF graphic")
    prs.add_argument("--ct", type=str, help="Base CT file to plot")
    prs.add_argument("--fasta", type=str, help="Fasta sequence file")
    prs.add_argument("--refct", type=str, help="Reference CT to compare base ct against")
    prs.add_argument("--pairprob", type=str, help="Pairing prob. file in dotplot format. By default, arcs are drawn for the following probability intervals: [0.03,0.1], [0.1,0.3], [0.3,0.8], [0.8,1.0]. These intervals can be modified by passing thresholds as comma-separated values. For example, --pairprob file.dp,0.03,0.1,0.3,0.8,1.0 specifies the default intervals. At most 4 intervals can be plotted, but fewer intervals are allowed (e.g. --pairprob file.dp,0.1,0.3,0.8,1.0 will plot 3 intervals).")
    prs.add_argument("--corr", type=str, help="Correlation file in dotplot format. Default color thresholds are 0.02 and 0.06. Can be modified by passing thresholds as comma-separated values (e.g. --corr file.dp,0.01,0.07)")
    
    prs.add_argument("--ntshape", type=str, help="Color nucs by shape reactivty in provided shape/map file")  
    prs.add_argument("--profile", type=str, help="Plot reactivity profile on top from shape/map file")
    
    prs.add_argument("--dmsprofile",type=str, help='Normalize and plot DMS reactivity from profile file')

    prs.add_argument("--bottom", action='store_true', help="Plot arcs on the bottom")

    prs.add_argument("--title", type=str, default='', help='Figure title')
    prs.add_argument("--showGrid", action="store_true", help="plot a grid on axis ticks")
    prs.add_argument("--intDistance", type=str, help="Pass depth file for computing likely max interaction distance")
    prs.add_argument("--depthThresh", type=int,default=10000, help="Interaction depth threshold (Default=10,000)")
    prs.add_argument("--bound", type=str, help="comma separated bounds of region to plot (e.g. --bound 511,796)")

    args = prs.parse_args()


    # subparse the pairprob argument
    args.ppbins = None
    if args.pairprob:
        ppspl = args.pairprob.split(',')
        bins = None
        if len(ppspl) > 1:
            args.pairprob = ppspl[0]
            args.ppbins = map(float, ppspl[1:])


    # subparse the bounds argument
    if args.bound:
        args.bound = map(int, args.bound.split(','))
    
    # subparse the corr argument
    args.colorbins = None
    if args.corr:
        corrspl = args.corr.split(',')
        if len(corrspl) == 3:
            args.corr = corrspl[0]
            args.colorbins = [-1, -float(corrspl[2]), -float(corrspl[1]), 
                               0, float(corrspl[1]), float(corrspl[2]), 1]
        elif len(corrspl) != 1:
            exit("Invalid --corr argument :: {0}".format(arg.corr))

    return args





if __name__=="__main__":
 
    arg = parseArgs()
    
    if arg.refct and not arg.ct:
        exit("--refct is invalid without --ct")
    
    
    aplot = ArcPlot(title = arg.title, fasta=arg.fasta)

    panel = 1
    if arg.bottom:
        panel = -1

    if arg.ct:
        if arg.refct:
            aplot.compareCTs( RNAtools.CT(arg.refct), RNAtools.CT(arg.ct), panel=panel)
        else:
            aplot.addCT( RNAtools.CT(arg.ct, filterNC=False), panel=panel)
        panel *= -1


    if arg.pairprob:
        aplot.addDotPlot( RNAtools.DotPlot(arg.pairprob), scheme = 'pairprob', panel=panel, bins=arg.ppbins)
        panel *= -1


    if arg.corr:
        aplot.addDotPlot( RNAtools.DotPlot(arg.corr), scheme = 'corr', panel=panel, bins=arg.colorbins)


    if arg.intDistance:
        aplot.addInteractionDistance(arg.intDistance, arg.depthThresh, panel)
            

    if arg.ntshape:
        aplot.colorSeqByMAP(arg.ntshape)


    if arg.profile:  
        aplot.readProfile(arg.profile)
    
    if arg.dmsprofile:
        aplot.readDMSProfile(arg.dmsprofile)

    if arg.showGrid:
        aplot.grid = True

   
    aplot.writePlot( arg.outputPDF, bounds = arg.bound)








