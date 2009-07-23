import numpy as np
import Image
from rebin import congrid
from glob import glob


def calculate_global_offsets(fname):
    '''calculates global offsets for each stack from relative offsets
    to neighbouring stacks

    Parameters:
      fname : filename of offsets file
    '''
    
    f = open(fname)
    lines = f.readlines()
    unglobd = []
    globd = {}
    for line in lines:
        if not line[0] == '#' or line.strip() == '':
            bits = line.split()
            entry = { 'rel':bits[1],
                     'ref':bits[0] }
            unglobd.append(entry)
            if bits[2].isalpha():
                entry['rel'] += '|' + bits[2]
                entry['slice'] = bits[3]
                entry['offset'] = np.asarray([int(x) for x in bits[4:7]])
                # need to handle split stacks
            else:
                # next three entries are x y z
                entry['offset'] = np.asarray([int(x) for x in bits[2:5]])
    while len(unglobd) > 0:
        entry = unglobd.pop(0) # take first item in list
        if globd.has_key(entry['ref']):
            globd[entry['rel']] = globd[entry['ref']] + entry['offset']
        elif entry['ref'] == 'glob':
            globd[entry['rel']] = np.asarray((0,0,0))
        else:
            unglobd.append(entry)
    return globd 

def get_max_projections(gdir):
    '''returns paths to all max projections - files two directories deep,
    with name ending in _max.png'''
    stacks = glob(gdir+'*/*_max.png')
    return stacks

def generate_max_projections_montage(gdir, scale=1):
    '''constructs a maximum projection montage from maximum projections
    and offsets file'''
    offsets_file = gdir + 'stack_offsets.txt'
    global_offsets = calculate_global_offsets(offsets_file)
    stacks = get_max_projections(gdir)
    
    # match stack name with offset
    for short_code in global_offsets.keys():
        short_name = short_code.split('|')[0]
        i = 0
        found = False
        while not found:
            if stacks[i].find(short_name) <> -1:
                found = True
                global_offsets[short_code] = \
                                           [global_offsets[short_code], \
                                            stacks[i]]
            i += 1

    # calculate biggest size of array
    # assume 1024x1024 + biggest offset

    # find biggest x and biggest y offsets
    maxx, maxy = 0, 0
    minx = miny = 0
    for info in global_offsets.itervalues():
        maxx = np.maximum(maxx, info[0][0])
        maxy = np.maximum(maxy, info[0][1])
        minx = np.minimum(minx, info[0][0])
        miny = np.minimum(miny, info[0][1])

    bigsz = np.asarray((maxx + 1024 + np.abs(minx),
                        np.abs(miny) + 1024 + maxy)) / float( scale )
    
    big = np.zeros(np.hstack((bigsz,2)))
    for info in global_offsets.itervalues():
        #sts = ['upa', 'som', 'rta', 'upcd', 'rtb', 'rtc', 'dnr', 'dnba', \
        #       'dnc', 'dnl', 'dnga', 'dngb']
        #for info in [global_offsets[x] for x in sts]:
        dat = np.asarray(Image.open(info[1]))
        datsz = np.asarray(dat.shape)
        datnewshape = np.round(datsz / float( scale )).astype(int)
        dat = congrid(dat, datnewshape)
        startx = (-minx + info[0][0]) / float( scale )
        starty = (-miny + info[0][1]) / float( scale )
        
        big[startx:startx + datnewshape[0], \
            starty:starty + datnewshape[1],1] = dat
        big[:,:,0] = big.max(2)
    return big.max(2)
