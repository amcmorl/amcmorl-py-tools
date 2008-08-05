#####################
## Version History ##
#####################

clusterplotversion = '0.1'
clusterplotmodified = '08 September, 2006'
clusterplotmodifier = 'R. Padraic Springuel'
#Initial creation date.

clusterplotversion = '0.2'
clusterplotmodified = '09 September, 2006'
clusterplotmodifier = 'R. Padraic Springuel'
#Debugging of treebuild.  Wrote poptree.  Modified datasort to include the
#option to favor one side for the heavier branch.

clusterplotversion = '0.3'
clusterplotmodified = '10 September, 2006'
clusterplotmodifier = 'R. Padraic Springuel'
#Added datalabels and wholetree.  Added kwargs to buildtree.  Beta release.

clusterplotversion = '0.4'
clusterplotmodified = '10 September, 2006'
clusterplotmodifier = 'R. Padraic Springuel'
#Syntax fixes in datalabels and wholetree.

clusterplotversion = '0.5'
clusterplotmodifed = '11 September, 2006'
clusterplotmodifier = 'R. Padraic Springuel'
#Added GNU General Public License information.

clusterplotversion = '0.6'
clusterplotmidified = '13 September, 2006'
clusterplotmodifier = 'R. Padraic Springuel'
#Changed lisense to BSD-style.  Given that this package is based on Pycluster
#and matplotlib, both of which uses licenses which are more closely aligned with
#the BSd-style than the GPL this makes logical sense.  Also, if others find this
#prackage useful enough to warrant its inclusion in matplotlib, this license
#change allows that to happen.

clusterplotversion = '0.7'
clusterplotmodified = '15 September, 2006'
clusterplotmodifier = 'R. Padraic Springuel'
#Modified coords to allow for the use of alternate metrics to determine the
#y-coordinates.

clusterplotversion = '0.8'
clusterplotmodified = '15 September, 2006'
clusterplotmodifier = 'R. Padraic Springuel'
#Added partial directional control to treebuild.  The use of linestyle='steps'
#causes the horizontal dendrogram to look wrong, however.

clusterplotversion = '1.0'
clusterplotmodified = '16 September, 2006'
clusterplotmodifier = 'R. Padraic Springuel'
#Completed directional control.  treebuild no longer uses linestyle='steps'
#which opens up more formating options and fixes the problem with the look of
#the horizontal dendrogram.  Inversion of the diagram is now also possible.
#First Full Release.

clusterplotversion = '1.1'
clusterplotmodified = '20 September, 2006'
clusterplotmodifier = 'R. Padraic Springuel'
#Bug fix: coords was stuck in integer math for its horizontal coordinates
#calculations.  Changed to floating point math.

clusterplotversion = '1.2'
clusterplotmodified = '21 September, 2006'
clusterplotmodifier = 'R. Padraic Springuel'
#New wrapper function: fadetree, variant on poptree that fades the branches
#according to population rather than printing the numbers.

clusterplotversion = '1.3'
clusterplotmodified = '26 September, 2006'
clusterplotmodifier = 'Michael Sorich & R. Padraic Springuel'
#Bug fix: datalabels was raising an exception and didn't contain orientation
#control.

clusterplotversion = '1.3.1'
clusterplotmodified = '28 September, 2006'
clusterplotmodifier = 'R. Padraic Springuel'
#Syntax fix: two if statements had '=' instead of '=='

#####################
## Program Summary ##
#####################

#Visualization of clustering solutions, in particular dendrograms for
#hierarchical clustering.  These routines assume Pycluster was used to generate
#the clustering solutions and that matplotlib is being uses matplotlib to create
#the graphs.

#At the moment the dendrograms are constructed vertically with the data on the
#bottom.  It is my intention to add at some point routines that allow the
#orientation of the dendrogram to be in any cardinal direction.

#If these routines are used with interactive mode on, your program will slow
#down signifigantly for large dendrograms because this program relies on
#multiple calls of different matplotlib commands.  This package will work better
#if you use it with interactive mode off and issue a show() command after you've
#called all the routines from this package that you need.

#NOTE: The comment documentation will alwasy describe the dendrogram as if it
#was oriented vertically.  It is possible, however, to orient the dendrogram
#horizontally.  If you choose this option, then what the documentation describes
#as happening vertically will happen horizontally and vis versa.

##########################
## Liscense Information ##
##########################

###############################################################################
#Copyright (c) 2006, R. Padraic Springuel                                     #
#All rights reserved.                                                         #
#                                                                             #
#Redistribution and use in source and binary forms, with or without           #
#modification, are permitted provided that the following conditions are met:  #
#                                                                             #
#    * Redistributions of source code must retain the above copyright notice, #
#      this list of conditions and the following disclaimer.                  #
#    * Redistributions in binary form must reproduce the above copyright      #
#      notice, this list of conditions and the following disclaimer in the    #
#      documentation and/or other materials provided with the distribution.   #
#    * Neither the name of the University of Maine, Orono nor the names of    #
#      its contributors may be used to endorse or promote products derived    #
#      from this software without specific prior written permission.          #
#                                                                             #
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"  #
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    #
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE   #
#ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE     #
#LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR          #
#CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF         #
#SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS     #
#INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      #
#CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)      #
#ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE   #
#POSSIBILITY OF SUCH DAMAGE.                                                  #
###############################################################################


##########################
## Program Dependancies ##
##########################

#These routines require the following packages to work properly:
#   matplotlib
#   numpy
#   Pycluster (which is dependant on Numeric)

#######################################################################

from pylab import *
from Pycluster import *

#This routine sorts the data based on the information in the cluster tree so
#that when the tree is created there are no branch crossings.  The return is
#a list of the indecies in an order that will result in no branch crossings.
#The heavy option can be set to 'left' or 'right' and when set the larger group
#will always be placed on that side.  If heavy is set, then you can also specify
#a weight for each data point.
def datasort(tree,heavy=None,weight=None):
    if heavy != None:
        pop = clusterpop(tree,weight)
    if heavy == 'left':
        if pop[tree[-1].left] > pop[tree[-1].right]:
            order = [tree[-1].left,tree[-1].right]
        else:
            order = [tree[-1].right,tree[-1].left]
    elif heavy == 'right':
        if pop[tree[-1].left] > pop[tree[-1].right]:
            order = [tree[-1].right,tree[-1].left]
        else:
            order = [tree[-1].left,tree[-1].right]
    else:
        order = [tree[-1].left,tree[-1].right]
    for i in range(-(len(tree)),0):
        for j in range(len(order)):
            if i == order[j]:
                if heavy == 'left':
                    if pop[tree[abs(i)-1].left] > pop[tree[abs(i)-1].right]:
                        order[j] = tree[abs(i)-1].right
                        order.insert(j,tree[abs(i)-1].left)
                    else:
                        order[j] = tree[abs(i)-1].left
                        order.insert(j,tree[abs(i)-1].right)
                elif heavy == 'right':
                    if pop[tree[abs(i)-1].left] > pop[tree[abs(i)-1].right]:
                        order[j] = tree[abs(i)-1].left
                        order.insert(j,tree[abs(i)-1].right)
                    else:
                        order[j] = tree[abs(i)-1].right
                        order.insert(j,tree[abs(i)-1].left)
                else:
                    order[j] = tree[abs(i)-1].right
                    order.insert(j,tree[abs(i)-1].left)
    return order


#This routine calculates the population of each cluster in a hierarchical
#clustering solution and returns them as a list of populations.  Taking
#advantage of the fact that Pycluster uses negative numbers to identify clusters
#and positive numbers to identify data points, the first half of the array
#contains the individual data points while the second half of the array are the
#clusters listed from the last formed to first.  Thus pop[-1] (the last entry
#in the list) is the population of the cluster that Pycluster identifies as -1
#in the tree (i.e. the first cluster it made), pop[-2] (the second to last
#entry) is the population of the cluster that Pycluster identifies as -2 in
#the tree, etc (i.e. the second cluster it made).  As in Pycluster, weights
#are taken to be the number of times a particular data point occurs within
#the data (i.e. the population of that data point) but need not be whole
#numbers.  If no weights are given, then the function assumes a weight of 1 for
#all data points.
def clusterpop(tree,weight=None):
    from numpy import zeros
    pop = zeros(len(tree)*2+1)
    if weight != None:
        if len(weight) == len(tree)+1:
            for i in range(len(weight)):
                pop[i] = weight[i]
        else:
            print 'There must be the same number of weights as there are data points.'
            return
    else:
        for i in range(len(tree)+1):
            pop[i] = 1
    for i in range(len(tree)):
        pop[-1-i] = pop[tree[i].left] + pop[tree[i].right]
    return pop


#The first step in plotting a dendrogram, this routine calculates the
#coordinates of each data point and cluster.  As with clusterpop, the fact
#that Pycluster uses negative numbers to identify its clusters is exploited
#to identify the coordinate tuples.  Thus coords[-1] is the coordinates of the
#first cluster created, coords[-2] is the coordinates of the second cluster
#created, etc.
#By default the data are placed at zero height.  This can be changed by setting
#the zero argument to a different value.
#Under default circumstances, the vertical coordinates for a cluster are set
#by the distance between its two children.  However, this behavior can be
#changed by providing a distalt arguement, which should be a list of numbers
#the same length as the tree.  This option also recognizes 'order' as a valid
#values for distalt.  If 'order' is specified, then the vertical coordinate is
#cluster's ordinal position in the tree (i.e. the first cluster made gets a
#height of 1, the second 2, etc.).
#Under default circumstances, the horizontal coordinates for a cluster are a
#weighted average of its child clusters with the weights being determined by
#the population of those child clusters.  This behavior can be changed so that
#the horizonzal coordinates for a cluster are a simple average of its child
#clusters by setting the sym argument to True.  Note, however, that if more than
#two clusters get tied together into one cluster at the same cluster distance
#this joining will not look symetric since Pycluster has to make this joining
#in several steps.
def coordinates(tree,zero=0,distalt=None,heavy=None,weight=None,sym=False):
    from numpy import zeros
    coords = list(zeros(len(tree)*2+1))
    a = datasort(tree,heavy,weight)
    if distalt == 'even':
        split = splits(tree)
    for i in range(len(a)):
        coords[a[i]] = [float(i),zero]
    if not sym:
        pop = clusterpop(tree,weight)
    for i in range(len(tree)):
        if distalt == None:
            y = tree[i].distance
        elif distalt == 'order':
            y = i+1
        else:
            y = distalt[i]
        if not sym:
            x = (pop[tree[i].right]*coords[tree[i].right][0]+pop[tree[i].left]*coords[tree[i].left][0])/pop[-1-i]
        else:
            x = (coords[tree[i].right][0]+coords[tree[i].left][0])/2
        coords[-1-i] = [x,y]
    return coords


#The workhorse function of this program, this routine takes the coordinates and
#the tree and uses them to construct the actual dendrogram.
#The mask argument allows you to hide clusters and data points.  It should be
#a 1-d list or array of length len(tree)*2+1 containing 0's  (do not plot) and
#1's (plot).  The first len(tree)+1 elements are the data points, the remainder
#are the clusters in reverse order (so mask[-1] corresponds to the first cluster
#made in the tree, which Pycluster identifies as -1).  Both ends of a conection
#must be unmasked in order for it to be drawn.
#By calling this function multiple times with different mask and line settings
#it is possible to set different branches to different colors and line settings.
#The last argument, p, represents the portion of the vertical plot space that
#the tree will take up.  The remainder is split 3:1 above and below the tree.
#Setting p to be something below 1 ensures that the entire vertical extent of
#tree will be visible.  The horizontal extents of the tree are handled by
#setting the x-axis limits to run from -1 to len(tree)+2.  This extends the
#x-axis one point beyond the data at either end.
#If invert is set to False (the default setting) then the vertical axis will be
#plotted in the traditional manner with small numbers on the bottom and large
#ones on top.  If invert is set to True then the vertical axis will be plotted
#in reverse, with large numbers on the bottom and small ones on top.  Note that
#invert controls the orientation of the axes directly, not the orientation of
#the tree.  As a result, if you chose a distalt in coords that would naturally
#reverse the direction of the tree, setting invert to True will set the tree
#direction back to normal (a single cluster on the top).
def treebuild(coords,tree,mask=None,orient='v',invert=False,line='b-',p=0.95,**kwargs):
    from numpy import array
    if mask == None:
        mask = ones(len(coords))
    hold(True)
    for i in range(len(tree)):
        if mask[-1-i] and mask[tree[i].left]:
            if orient == 'v':
                x = [coords[-1-i][0],coords[tree[i].left][0],coords[tree[i].left][0]]
                y = [coords[-1-i][1],coords[-1-i][1],coords[tree[i].left][1]]
            elif orient == 'h':
                y = [coords[-1-i][0],coords[tree[i].left][0],coords[tree[i].left][0]]
                x = [coords[-1-i][1],coords[-1-i][1],coords[tree[i].left][1]]
            else:
                print 'Invalid orientation.'
                return
            plot(x,y,line,**kwargs)
        if mask[-1-i] and mask[tree[i].right]:
            if orient == 'v':
                x = [coords[-1-i][0],coords[tree[i].right][0],coords[tree[i].right][0]]
                y = [coords[-1-i][1],coords[-1-i][1],coords[tree[i].right][1]]
            elif orient == 'h':
                y = [coords[-1-i][0],coords[tree[i].right][0],coords[tree[i].right][0]]
                x = [coords[-1-i][1],coords[-1-i][1],coords[tree[i].right][1]]
            plot(x,y,line,**kwargs)
    if orient == 'v':
        xlim(-1,len(tree)+2)
        ymin = coords[0][1]
        ymax = coords[0][1]
        for i in range(len(coords)):
            if coords[i][1] < ymin:
                ymin = coords[i][1]
            elif coords[i][1] > ymax:
                ymax = coords[i][1]
        if invert:
            ylim(ymax+(ymax-ymin)*(1.-p)*3/4,ymin-((ymax-ymin)*(1.-p)/4))
        else:
            ylim(ymin-((ymax-ymin)*(1.-p)/4),ymax+(ymax-ymin)*(1.-p)*3/4)
    elif orient == 'h':
        ylim(-1,len(tree)+2)
        xmin = coords[0][1]
        xmax = coords[0][1]
        for i in range(len(coords)):
            if coords[i][1] < xmin:
                xmin = coords[i][1]
            elif coords[i][1] > xmax:
                xmax = coords[i][1]
        if invert:
            xlim(xmax+(xmax-xmin)*(1.-p)*3/4,xmin-((xmax-xmin)*(1.-p)/4))
        else:
            xlim(xmin-((xmax-xmin)*(1.-p)/4),xmax+(xmax-xmin)*(1.-p)*3/4)
    return


#This function labels each cluster with the specified string.  The use of
#mask can prevent some clusters from being labeled.  The text options can be
#specified via fontdict and/or keyword arguments as normal for matplotlib.
def clusterlabels(coords,labels,mask=None,fontdict=None,**kwargs):
    from numpy import ones
    if mask == None:
        mask = ones(len(coords))
    for i in range(len(coords)):
        if mask[i] == 1:
            text(coords[i][0],coords[i][1],labels[i],fontdict=fontdict,**kwargs)
    return


#Thus function uses the x-axis tick marks to label the data points.  Of course
#they could also be labeled using the clusterlabels function because each data
#point can also be considered a cluster, but this makes the placement of data
#labels seperately from cluster labels easier, especially if you only want to
#label the data and not the clusters.
def datalabels(tree,dlabels,heavy=None,weight=None,mask=None,orient='v',fontdict=None,**kwargs):
    from numpy import ones
    if mask == None:
        mask = ones(len(dlabels))
    a = datasort(tree,heavy,weight)
    dl = []
    for i in range(len(dlabels)):
        if mask[i]:
            dl.append(dlabels[a[i]])
        else:
            dl.append('')
    if orient == 'v':
        xticks(range(len(dlabels)),dl,fontdict=fontdict,**kwargs)
    if orient == 'h':
        yticks(range(len(dlabels)),dl,fontdict=fontdict,**kwargs)
    return


#######################
## Wrapper Functions ##
#######################

#The wrapper functions below allow the quick creation of dendrograms.  They do
#that, however, by limiting the options avaliable to manipulate the appearence
#of the dendrogram.  For full control over the appearence of the dendrogram,
#the above functions should be used independantly.

#To minimize conflicts, the color of the lines in the tree cannot be specified
#with a matplotlib keyword argument but must instead be specified with the
#specific color argument.  The keyword arguments are reserved for text
#text properties for the labels.


#This wrapper function simply creates the whole tree and places the data labels
#if they are provided.
def wholetree(tree,dlabels=None,heavy=None,weight=None,line='b',sym=False,p=0.95,fontdict=None,**kwargs):
    coords = coordinates(tree,0,None,heavy,weight,sym)
    treebuild(coords,tree,None,'v',False,line,p)
    if dlabels != None:
        datalabels(tree,dlabels,heavy,weight,'v',fontdict=fontdict,**kwargs)
    return

    
#This wrapper function simplifies the creatation of a dendrogram where each
#cluster is labeled by its population.  The optional lowerlimit argument sets
#the maximum size of the clusters that are plotted.  Identical data points
#(i.e. those with 0 joining distance) are always hidden.
def poptree(tree,heavy=None,weight=None,line='b',sym=False,p=0.95,lowerlimit=1,fontdict=None,**kwargs):
    from numpy import ones
    pop = clusterpop(tree)
    mask = ones(len(pop))
    labels = []
    for i in range(len(pop)):
        if pop[i] < lowerlimit:
            mask[i] = 0
        labels.append(str(pop[i]))
    for i in range(len(tree)):
        if tree[i].distance == 0:
            mask[tree[i].right] = 0
            mask[tree[i].left] = 0
    coords = coordinates(tree,0,None,heavy,weight,sym)
    treebuild(coords,tree,mask,'v',False,line,p)
    clusterlabels(coords,labels,mask,fontdict=fontdict,**kwargs)
    return

#A variation on poptree that uses alpha transparency to fade branches according
#to population instead of printing the numbers on the tree.
def fadetree(tree,heavy=None,weight=None,line='b',sym=False,p=0.95,fontdict=None,**kwargs):
    from numpy import zeros
    pop = clusterpop(tree)
    coords = coordinates(tree,0,None,heavy,weight,sym)
    mask1 = zeros(len(pop))
    mask2 = zeros(len(pop))
    mask3 = zeros(len(pop))
    mask4 = zeros(len(pop))
    mask5 = zeros(len(pop))
    mask6 = zeros(len(pop))
    mask7 = zeros(len(pop))
    mask8 = zeros(len(pop))
    mask9 = zeros(len(pop))
    mask0 = zeros(len(pop))
    for i in range(len(pop)):
        if float(pop[i])/max(pop) < 0.1:
            mask1[i] = 1
        elif float(pop[i])/max(pop) < 0.2:
            mask2[i] = 1
        elif float(pop[i])/max(pop) < 0.3:
            mask3[i] = 1
        elif float(pop[i])/max(pop) < 0.4:
            mask4[i] = 1
        elif float(pop[i])/max(pop) < 0.5:
            mask5[i] = 1
        elif float(pop[i])/max(pop) < 0.6:
            mask6[i] = 1
        elif float(pop[i])/max(pop) < 0.7:
            mask7[i] = 1
        elif float(pop[i])/max(pop) < 0.8:
            mask8[i] = 1
        elif float(pop[i])/max(pop) < 0.9:
            mask9[i] = 1
        else:
            mask0[i] = 1
    for i in range(len(tree)):
        if mask1[tree[i].left] or mask1[tree[i].right]:
            mask1[-1-i] = 1
        if mask2[tree[i].left] or mask2[tree[i].right]:
            mask2[-1-i] = 1
        if mask3[tree[i].left] or mask3[tree[i].right]:
            mask3[-1-i] = 1
        if mask4[tree[i].left] or mask4[tree[i].right]:
            mask4[-1-i] = 1
        if mask5[tree[i].left] or mask5[tree[i].right]:
            mask5[-1-i] = 1
        if mask6[tree[i].left] or mask6[tree[i].right]:
            mask6[-1-i] = 1
        if mask7[tree[i].left] or mask7[tree[i].right]:
            mask7[-1-i] = 1
        if mask8[tree[i].left] or mask8[tree[i].right]:
            mask8[-1-i] = 1
        if mask9[tree[i].left] or mask9[tree[i].right]:
            mask9[-1-i] = 1
        if mask0[tree[i].left] or mask0[tree[i].right]:
            mask0[-1-i] = 1
    treebuild(coords,tree,mask1,'v',False,line,p,alpha=0.1)
    treebuild(coords,tree,mask2,'v',False,line,p,alpha=0.2)
    treebuild(coords,tree,mask3,'v',False,line,p,alpha=0.3)
    treebuild(coords,tree,mask4,'v',False,line,p,alpha=0.4)
    treebuild(coords,tree,mask5,'v',False,line,p,alpha=0.5)
    treebuild(coords,tree,mask6,'v',False,line,p,alpha=0.6)
    treebuild(coords,tree,mask7,'v',False,line,p,alpha=0.7)
    treebuild(coords,tree,mask8,'v',False,line,p,alpha=0.8)
    treebuild(coords,tree,mask9,'v',False,line,p,alpha=0.9)
    treebuild(coords,tree,mask0,'v',False,line,p,alpha=1.0)
    return
    

def ClusterPlotversion():
    print 'ClusterPlot.py'
    print 'Created 08 September, 2006'
    print 'by R. Padraic Springuel'
    print 'Version %s modified %s' % (clusterplotversionnum,clusterplotmodified)
    print 'Most recent modification by %s' % clusterplotmodifier
    return
