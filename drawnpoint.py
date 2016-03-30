""" draw n-point functions in matplotlib
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.collections import PatchCollection

import numpy as np

class Canvas(object):
    def __init__(self, xsize=5, ysize=4):
        self.plt = plt.figure(figsize=(xsize, ysize))
        self.ax = plt.axes()
        self.arstyle = "-|>"
        self.arprops = dict(arrowstyle=self.arstyle, connectionstyle="arc3")
        self.patches = []

    def draw_triangle(self, xv1, xv2, xv0=np.array([0.,0.0])):
        """
        draw a triangle with arrows;
        each argument is a vector i.e. xv1 = (x1, y1) etc.
        The third side then, by definition is
        xv3 = -xv1 - xv2 = (-x1-x2, -y1-y2)
        xv0 defines the origin; default is (0.0, 0.0)
        """
        # draw arrows
        xvalues = np.array([xv0[0], xv1[0], xv2[0]])
        yvalues = np.array([xv0[1], xv1[1], xv2[1]])
        minx, maxx, miny, maxy = min(xvalues), max(xvalues), min(yvalues), max(yvalues)
        ds = np.sqrt((maxx-minx)**2.0 + (maxy-miny)**2.0)

        arstyle = self.arstyle
        arprops = dict(arrowstyle=self.arstyle, connectionstyle="arc3")

        self.ax.annotate("", xytext=(xv0[0], xv0[1]), xycoords='data', xy=(xv1[0], xv1[1]), textcoords='data', arrowprops=self.arprops)
        self.ax.annotate("", xytext=(xv1[0], xv1[1]), xycoords='data', xy=(xv2[0], xv2[1]), textcoords='data', arrowprops=self.arprops)
        self.ax.annotate("", xytext=(xv2[0], xv2[1]), xycoords='data', xy=(xv0[0], xv0[1]), textcoords='data', arrowprops=self.arprops)

        plt.xlim(minx-ds/10., maxx+ds/10.)
        plt.ylim(miny-ds/10., maxy+ds/10.)

    def draw_quad(self, xv1, xv2, xv3, xv0=np.array([0., 0.])):
        """
        draw a quadrilateral with arrows;
        each argument is a vector i.e. xv1 = (x1, y1) etc.
        The third side then, by definition is
        xv3 = -xv1 - xv2 -xv3 = (-x1-x2-x3, -y1-y2-y3)
        xv0 defines the origin; default is (0.0, 0.0)
        """
        # draw arrows
        xvalues = np.array([xv0[0], xv1[0], xv2[0], xv3[0]])
        yvalues = np.array([xv0[1], xv1[1], xv2[1], xv3[1]])
        minx, maxx, miny, maxy = min(xvalues), max(xvalues), min(yvalues), max(yvalues)
        ds = np.sqrt((maxx-minx)**2.0 + (maxy-miny)**2.0)

        self.ax.annotate("", xytext=(xv0[0], xv0[1]), xycoords='data', xy=(xv1[0], xv1[1]), textcoords='data', arrowprops=self.arprops)
        self.ax.annotate("", xytext=(xv1[0], xv1[1]), xycoords='data', xy=(xv2[0], xv2[1]), textcoords='data', arrowprops=self.arprops)
        self.ax.annotate("", xytext=(xv2[0], xv2[1]), xycoords='data', xy=(xv3[0], xv3[1]), textcoords='data', arrowprops=self.arprops)
        self.ax.annotate("", xytext=(xv3[0], xv3[1]), xycoords='data', xy=(xv0[0], xv0[1]), textcoords='data', arrowprops=self.arprops)

        plt.xlim(minx-ds/10., maxx+ds/10.)
        plt.ylim(miny-ds/10., maxy+ds/10.)

    def draw_square(self, xv0=np.array([0., 0.]), x=.5, ec="b"):
        """
        draw a square starting at the given position
        """
        sq = mpatches.Rectangle(xv0, width=x, height=x, fc='white', ec=ec, fill=False)
        self.ax.add_artist(sq)
        return sq

    def draw_squeezed_four_point(self):
        self.draw_quad([1.5,1.5], [0.5,4], [-0.1,0.4])
        self.ax.text(1, 3, r"$\vec{k}_2$")
        self.ax.text(0.7, 1.1, r"$\vec{k}_1$")
        self.ax.text(0.0, 2.2, r"$\vec{k}_3$")
        self.ax.text(0.05, 0.4, r"$\vec{q}$")
        self.ax.axis('off')

    def draw_squeezed_three_point(self):
        self.draw_triangle([1.5,1.5], [-0.1, 0.4])
        self.ax.text(0.7, 0.5, r"$-\vec{k}-\vec{q}$")
        self.ax.text(0.7, 1.1, r"$\vec{k}$")
        self.ax.text(0.03, 0.2, r"$\vec{q}$")
        self.ax.axis('off')

    def draw_subvolumes(self, Lbox=16., Nsub=4, ec="black"):
        dx = Lbox/Nsub
        for i in range(Nsub):
            stx = i*dx
            for j in range(Nsub):
                sty = j*dx
                self.draw_square(xv0=np.array([stx, sty]), x=dx, ec=ec)

        plt.xlim(-1, Lbox+1)
        plt.ylim(-1, Lbox+1)
