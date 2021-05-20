import sys
import math as m
import numpy as np
from PyQt4 import QtGui, QtCore
from GuiMainWindow import GuiMainWindow
from GuiTrailExecutor import GuiTrailExecutor

#-----------------------------------------------------------------------------
def main():
    app = QtGui.QApplication(sys.argv)
    the_main_wnd = GuiMainWindow()

    #trail_executor = GuiTrailExecutor(the_main_wnd)
    #trail_thread = QtCore.QThread()
    #trail_executor.moveToThread(trail_thread)
    #trail_thread.started.connect(trail_executor.run)
    #trail_thread.start()

    #workerThread = QtCore.QThread()
    #trail_executor.moveToThread(workerThread)
    #workerThread.start()

#    trail_executor = GuiTrailExecutor(the_main_wnd)
#    trail_executor.start()

    the_main_wnd.show()
    app.exec_()

#-----------------------------------------------------------------------------
if __name__ == '__main__':
    main()

#=============================== END OF FILE =================================
