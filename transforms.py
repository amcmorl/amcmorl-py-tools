import numpy


def rotate90degs(arr, number=1, direction='cw'):
    '''Uses transpose and index flipping if necessary
    to rotate arrays 90 degs in given axis pair

    currently only rotating in axes 0 and 1 is supported

    direction can be given as c, 0 or 1 for clock-wise
    and -1 or ccw for counter-clock-wise'''

    r = numpy.rank(arr)
    trlist = [1,0] + range(2,r)

    for i in xrange(number):
        # base 90 deg CW rotation
        arr = arr.transpose(trlist)
        if direction in ['cw', 0, 1]:
            arr = arr[:,::-1,...]
        elif direction in ['ccw',-1]:
            arr = arr[::-1,...]
        else:
            raise ValueError("unknown direction specified")

    return arr


def flip(arr, direction='h'):
    '''Uses index flipping to flip arrays

    direction can be either h or 0 for horizontal (flipping in axis 0)
    or v or 1 for vertical (flipping in axis 1)'''

    if direction in ['h', 0]:
        return arr[::-1,]
    elif direction in ['v', 1]:
        return arr[...,::-1]
    else:
        raise ValueError("unknown direction specified")
