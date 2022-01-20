import pyqtgraph as pg
import numpy as np
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtCore import QSize, Qt
import pyqtgraph.opengl as gl

xVals=(1,2,3,4)
yVals=(5,6,7,8)
zVals=(9,10,11,12)

pos = np.array([xVals, yVals, zVals])


# Initialize Qt (once per application)
app = QtWidgets.QApplication([])
# Define top-level widget to hold everything
w = QtWidgets.QWidget()
w.setAttribute(Qt.WA_StyledBackground, True)
w.setStyleSheet('background-color: white;')
    # Create a plot and a text (widgets)
    #plotWidget = pg.PlotWidget()

# Make a widget for displaying 3D stuff
glvw = gl.GLViewWidget()
###glvw.setBackgroundColor('w')



# Create grids for the 3D object, rotate and scale
xgrid = gl.GLGridItem()
ygrid = gl.GLGridItem()
zgrid = gl.GLGridItem()
glvw.addItem(xgrid)
xgrid.setColor('k')
ygrid.setColor('k')
zgrid.setColor('k')
glvw.addItem(ygrid)
glvw.addItem(zgrid)
xgrid.rotate(90, 0, 1, 0)
ygrid.rotate(90, 1, 0, 0)
xgrid.scale(1, 1, 1)
ygrid.scale(1, 1, 1)
zgrid.scale(1, 1, 1)

# Populate the widget
colors = np.array([(1, 1, 1, 1) for i in range(len(xVals))])


p3d = gl.GLScatterPlotItem(pos=pos, color=(0.25,0.25,0.5,0.5), size=10)
p3d.setGLOptions('opaque')
glvw.addItem(p3d)

# Populate the plot
#plotWidget.plot(xVals, yVals, zVals, pen=3)

text = QtWidgets.QLineEdit('Textorama')
# Create a grid layout to manage the widgets size and position
layout = QtWidgets.QGridLayout()
# Set this layout to be applied to the top-level widget that holds everything
w.setLayout(layout)
# Add widgets to the grid
#layout.addWidget(text, 0, 0)
#layout.addWidget(plotWidget, 0, 1, 3, 1)
layout.addWidget(glvw, 0, 0)
layout.setColumnMinimumWidth(0,500)
layout.setRowMinimumHeight(0,500)

#glvw.sizeHint = lambda: QSize(500, 500)
# Display the top-level widget in a new window
glvw.setBackgroundColor('w')





w.show()

# Start the Qt event loop
app.exec_()

