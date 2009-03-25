import numpy as np
import matplotlib.pyplot as plt

dashes = [[20, 2],
          [2, 2],
          [20, 2, 2, 2],
          [20, 2, 20, 2],
          [20, 2, 2, 2, 2, 2],
          [20, 2, 20, 2, 2, 2],
          [8, 2, 2, 2, 2, 2],
          [8, 2, 2, 2]]

colours = [[0., 0., 0.],
           [0., 0., 0.],
           [0.5, 0.5, 0.5],
           [0.5, 0.5, 0.5]]

class LineCycler():
    def __init__(self):
        self.c = 0
        self.d = 0

    def __call__(self, what='dashes'):
        if what == 'dashes' or what == 'd':
            n_styles = len(dashes)
            style = dashes[self.d % n_styles]
            self.d += 1
        else:
            n_styles = len(colours)
            style = colours[self.c % n_styles]
            self.c += 1
        return style

class FigNumer():
    def __init__(self):
        self.next_num = 0

    def __call__(self):
        next_num = self.next_num
        self.next_num += 1
        return next_num
