import numpy as np
from PyQt5 import QtWidgets, QtCore, QtOpenGL
from pyqtgraph.Qt import QtCore, QtGui
from PyQt5.QtCore import QSize, Qt
import pyqtgraph.opengl as gl
import pyqtgraph as pg
import OpenGL.GL as ogl


from cartesian import *
x_cart, y_cart, z_cart, times = parse_cartesian(path='/run/timeshift/backup/IOCB/cartesian/cart.xvg')
x_com, y_com, z_com = parse_com(path='/run/timeshift/backup/IOCB/cartesian/com.xvg')

x_cart_2, y_cart_2, z_cart_2, times = parse_cartesian(path='/run/timeshift/backup/IOCB/cartesian/trp_pdz.xvg')
x_com_2, y_com_2, z_com_2 = parse_com(path='/run/timeshift/backup/IOCB/cartesian/trp_pdz_com.xvg')

vectors = get_vectors(x_cart, y_cart, z_cart, x_com, y_com, z_com)
vectors_2 = get_vectors(x_cart_2, y_cart_2, z_cart_2, x_com_2, y_com_2, z_com_2)

# https://stackoverflow.com/questions/56890547/how-to-add-axis-features-labels-ticks-values-to-a-3d-plot-with-glviewwidget


"""
vectors = {
    'atom 1' : [(1,2,3), (1.2, 2.3, 3.1), (4,5,3), (0.5, 0.6, 0.8)],
    'atom 2' : [(4,7,8), (1, 5.5, 6), (9.7, 8.7, 6.7), (8.2, 2.2, 2.3)]
}
"""

# Initialize Qt (once per application)
app = QtGui.QApplication([])
# Define top-level widget to hold everything
w = gl.GLViewWidget()
w.resize(1920, 1080)
w.setAttribute(Qt.WA_StyledBackground, True)
w.setStyleSheet('background-color: white;')
# Create a grid layout to manage the widgets size and position
layout = QtWidgets.QGridLayout()
# Set this layout to be applied to the top-level widget that holds everything
w.setLayout(layout)

max_columns = 6
row_iter = 0
col_iter = 0

for vector_key_1, vector_key_2 in zip(list(vectors.keys())[:50], list(vectors.keys())[:10]):
    if col_iter == max_columns:
        col_iter = 0
        row_iter += 1

    pos = np.array([[x,y,z] for x,y,z in vectors[vector_key_1]])
    pos_2 = np.array([[x, y, z] for x, y, z in vectors[vector_key_2]])

    glvw = gl.GLViewWidget()
    label = QtWidgets.QLabel()
    label.setText(vector_key_1)
    xgrid = gl.GLGridItem()
    ygrid = gl.GLGridItem()
    zgrid = gl.GLGridItem()

    xgrid.setColor('k')
    ygrid.setColor('k')
    zgrid.setColor('k')
    #glvw.addItem(xgrid)
    #glvw.addItem(ygrid)
    #glvw.addItem(zgrid)
    xgrid.rotate(90, 0, 1, 0)
    ygrid.rotate(90, 1, 0, 0)
    xgrid.scale(5, 5, 5)
    ygrid.scale(5, 5, 5)
    zgrid.scale(5, 5, 5)







    #p3d = gl.GLScatterPlotItem(pos=pos, color=(0.25, 0.25, 0.5, 0.5), size=10)

    p3d = gl.GLLinePlotItem(pos=pos, color=(0.25, 0.25, 0.5, 0.5))
    glvw.addItem(p3d)
    #p3d_2 = gl.GLScatterPlotItem(pos=pos_2, color=(0.75, 0.75, 0.75, 0.75))
    p3d_2 = gl.GLLinePlotItem(pos=pos_2, color=(0.10, 0.75, 0.75, 0.75))
    glvw.addItem(p3d_2)


    p3d.setGLOptions('opaque')
    p3d_2.setGLOptions('opaque')

    axis = gl.GLAxisItem(glOptions="opaque")
    axis.setGLOptions('opaque')
    axis.setSize(x=10, y=10, z=10)

    glvw.addItem(axis)


    text = gl.GLTextItem(text=vector_key_1, color='k')
    glvw.addItem(text)

    # Create axis labels
    xlabel = gl.GLTextItem(text='X', color='k', pos=(2, 0.1, 0.1))
    ylabel = gl.GLTextItem(text='Y', color='k', pos=(0.1, 2, 0.1))
    zlabel = gl.GLTextItem(text='Z', color='k', pos=(0.1, 0.1, 2))



    glvw.addItem(xlabel)
    glvw.addItem(ylabel)
    glvw.addItem(zlabel)
    glvw.setBackgroundColor('w')
    glvw.opts['distance'] = 5.0
    glvw.opts['fov'] = 120.0
    glvw.opts['elevation'] = 45
    glvw.opts['azimuth'] = 45
    ###
    layout.addWidget(glvw, row_iter, col_iter)
    layout.setColumnMinimumWidth(col_iter, 50)
    layout.setRowMinimumHeight(row_iter, 50)

    col_iter += 1

# Show window

w.show()
w.setWindowTitle('Traj')
# Start the Qt event loop
app.exec_()
