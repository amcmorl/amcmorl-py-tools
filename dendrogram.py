import numpy as n
import pyx
import neuron.path as np
from glob import glob

class Branch:
    def __init__(self, length, n_children=2):
        self.length = length
        self.children = []
        self.n_children = n_children

    def add_child(self, branch):
        self.children.append(branch)
        self.n_children -= 1

    def print_branch(self, order=0):
        if order > 0:
            for i in xrange(order):
                print "   ",
            down = "%d" % (order)
        else:
            down = ""
        print "%s--- l = %.1f ---" % (down, self.length), 
        if len(self.children) == 0:
            print "x"
        else:
            print ""
        order += 1
        for child in self.children:
            child.print_branch(order)

    def draw_branch(self, canvas, start_pt, scale):
        h_scale = 0.01
        sx, sy = start_pt
        fx = sx + self.length * h_scale
        spup = sy - scale
        spdown = sy + scale

        line = pyx.path.line(sx, sy, fx, sy)
        canvas.stroke(line)
            
        if len(self.children) > 0:
            newscale = 0.45 * scale
            spacer = pyx.path.line(fx, spup, fx, spdown)
            canvas.stroke(spacer)
            self.children[0].draw_branch(canvas, (fx, spup), newscale)
            self.children[1].draw_branch(canvas, (fx, spdown), newscale)


def draw_tree(tree):
    primespace = 4.
    rho = 0.48

    nprimes = len(tree.children)
    c = pyx.canvas.canvas()

    # draw soma line (vertical)
    somaline = pyx.path.line(0, 0, 0, (nprimes - 1) * primespace)
    c.stroke(somaline)
    for i, child in enumerate(tree.children):
        child.draw_branch(c, (0, i * primespace), \
                          primespace * rho * 0.5)
    return c

    
def vert_spacing(k, rho):
    loc = [1]
    for n in xrange(1, k):
        loc.append(loc[-1] - rho**n)
    return loc


def parse_points(pts, children=None, dists=None):
    '''uses a watershed (forwards only) parsing of the dendritic tree
    to capture all branches and associate correctly with branch points'''
    
    if children == None:
        children = np.construct_children(pts)   
    if dists == None:
        dists = np.cumulative_dist(pts, children)
    npts = pts[0].size
    parents = pts[-1]

    queue = [] # points/directions to check

    somapts = n.where(parents == -1)[0]
    for i in xrange(somapts.size):
        queue.append( [somapts[i], 0.] ) # seed soma pts

    soma = Branch(0., n_children=somapts.size)
    curbranch = [soma]

    # while points in queue
    while queue:
        # check first point in queue
        tpt, refdist  = queue.pop()

        if children[tpt] == -1:
            # reached end of dendrite
            bpt = Branch(dists[tpt] - refdist)
            curbranch[-1].add_child(bpt)
            if curbranch[-1].n_children == 0:
                del curbranch[-1]
        else:
            if children[tpt] > 0:
                # branch point
                # 1) deal with last branch
                bpt = Branch(dists[tpt] - refdist)
                curbranch[-1].add_child(bpt)
                if curbranch[-1].n_children == 0:
                    del curbranch[-1]

                # 2) deal with next branch
                curbranch.append(bpt)
                queue.append( [children[tpt], refdist] )
                refdist = dists[tpt]
                
            if parents[tpt + 1] == 0:
                queue.append( [tpt+1, refdist] )

    return soma

def do_all(all_dir):
    cells = glob(all_dir + '*/*/')
    for cell in cells:
        pts = np.read_points(cell + 'global_coords.txt')
        tree = parse_points(pts)
        c = draw_tree(tree)
        c.writePDFfile(cell + 'dendrogram.pdf')
        c.writeEPSfile(cell + 'dendrogram.eps')
