import sys
import math as m
import numpy as np
from PyQt4 import QtGui, QtCore
from CSubd import get_angle, double_polygon, circle_avg
from GuiMainWindow import GuiMainWindow
from GuiTrailExecutor import GuiTrailExecutor

#-----------------------------------------------------------------------------
def main():
    app = QtGui.QApplication(sys.argv)
    the_main_wnd = GuiMainWindow()
    trail_executor = GuiTrailExecutor(the_main_wnd)
    trail_executor.start()
    the_main_wnd.show()
    app.exec_()

#-----------------------------------------------------------------------------
if __name__ == '__main__':
    main()

#=============================== END OF FILE =================================
