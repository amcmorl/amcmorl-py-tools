import tables
import numpy
import os
from spine.Spineclass import stage_initd, stage_userptd, \
     stage_measured, stage_fitd
import spine
import pylab

defname = '/home/amcmorl/data/spines/spines_db.h5'

class SpineData(tables.IsDescription):
    name = tables.StringCol(16)
    source_cell = tables.StringCol(32)
    source_stack = tables.StringCol(16)
    source_offset = tables.IntCol(shape=3)
    image_size = tables.IntCol(shape=3)
    image = tables.FloatCol(shape=(30,30,30))
    psf_fwhm = tables.FloatCol(shape=3)
    pixel_spacing = tables.FloatCol(shape=3)
    living = tables.BoolCol()

    head = tables.IntCol(shape=3)
    neck = tables.IntCol(shape=3)
    origin = tables.IntCol(shape=3)
    dend2 = tables.IntCol(shape=3)
    head_int = tables.FloatCol()
    neck_int = tables.FloatCol()

    dend_diam = tables.FloatCol()
    neck_length = tables.FloatCol()
    neck_diam = tables.FloatCol()
    head_diam = tables.FloatCol()
    angle_toz = tables.FloatCol()
    angle_roundz = tables.FloatCol()

def tofit(arr):
    res = numpy.zeros((30,30,30), dtype=arr.dtype)
    res[0:arr.shape[0], 0:arr.shape[1], 0:arr.shape[2]] = arr
    return res


def create_db(fname):
    spinedb = tables.openFile(fname, mode='w', \
                              title='spines database')
    spine_table = spinedb.createTable(spinedb.root, 'spine_table', \
                                      SpineData, "data about spine")
    return spinedb


class spinedb:
    def __init__(self, fname=None):
        if fname == None:
            fname = defname
        if os.path.isfile(fname):
            self.dbfile = tables.openFile(fname, mode='a')
        else:
            self.dbfile = create_db(fname)


    def add_spine(self, spine, living=False):
        spine_row = self.dbfile.root.spine_table.row
        spine_row['name'] = spine.name
        spine_row['source_cell'] = spine.big_name.split('/')[0]
        if len(spine.big_name.split('/')) > 1:
            spine_row['source_stack'] = spine.big_name.split('/')[1]
        else:
            spine_row['source_stack'] = spine.big_name
        spine_row['source_offset'] = spine.big_offset
        spine_row['image_size'] = spine.source_data.shape
        spine_row['image'] = tofit( spine.source_data )
        spine_row['living'] = living
        
        spine_row['head'] = spine.user_pts[0]
        spine_row['neck'] = spine.user_pts[1]
        spine_row['origin'] = spine.user_pts[2]
        spine_row['dend2'] = spine.user_pts[3]
        spine_row['head_int'] = spine.params['head_int']
        spine_row['neck_int'] = spine.params['neck_int']
    
        spine_row['dend_diam'] = spine.params['dend_diam']
        spine_row['neck_length'] = spine.params['neck_length']
        spine_row['neck_diam'] = spine.params['neck_diam']
        spine_row['head_diam'] = spine.params['head_diam']
        spine_row['angle_toz'] = spine.params['angle_toz']
        spine_row['angle_roundz'] = spine.params['angle_roundz']
        spine_row.append()
        self.dbfile.root.spine_table.flush()


    def delete_spine(self):
        '''need to be able to delete incorrect spine entries (rows)'''
        


    def close(self):
        self.dbfile.close()


    def __del__(self):
        self.close()


    def list_spines(self):
        spine_table = self.dbfile.root.spine_table
        return [x['name'] for x in spine_table.iterrows()]


def populate_spinedb():
    ''' populates spinedb from scratch from existing spine ncdf files
    '''
    bdir = '/home/amcmorl/data/spines/'

    lives = \
          ['040220/spd/spd_spa.nc', \
           '040525/bat25a_spa.nc', \
           '040525/bat25a_spb.nc', \
           '040525/bat25a_spc.nc', \
           '040617/cellb/live/five_spa_exdcn.nc', \
           '040707/liv/baspa.nc', \
           '040707/liv/baspb.nc', \
           '040707/liv/baspc.nc', \
           '040707/liv/baspd.nc']

    fixeds = \
           ['050822-2ca/050924/spa.nc', \
            '050822-2ca/050924/spb.nc', \
            '050822-2ca/050924/spc.nc', \
            '050822-2ca/050924/spd.nc', \
            '050822-2ca/050924/spe.nc', \
            '050822-2ca/050924/spf.nc', \
            '050822-2ca/050924/spg.nc', \
            '050822-2ca/050924/sph.nc', \
            '021125/mrsp/0608_analysis/mrspa.nc', \
            '021125/mrsp/0608_analysis/mrspb.nc', \
            '021125/mrsp/0608_analysis/mrspc.nc', \
            '040624/040903/uca.nc', \
            '040624/040903/ucb.nc', \
            '040624/040903/ucc.nc', \
            '040624/040903/ucd.nc', \
            '040624/040903/uce.nc', \
            '040624/040903/ucf.nc']

    spdb = spinedb()
    for indi in lives:
        print "### Importing %s" % (bdir + indi)

        if os.path.isfile(bdir+indi):
            sp = spine.Spine(bdir + indi)
            spdb.add_spine(sp, living=True)
        else:
            print "Not a file!"

    for indi in fixeds:
        print "### Importing %s" % (bdir + indi)

        if os.path.isfile(bdir+indi):
            sp = spine.Spine(bdir + indi)
            spdb.add_spine(sp, living=False)
        else:
            print "Not a file!"
