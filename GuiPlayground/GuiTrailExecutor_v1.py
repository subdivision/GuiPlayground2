from PyQt4 import QtGui, QtCore
from GuiMainWindow import GuiMainWindow
import time

class GuiTrailExecutor(QtCore.QThread):
    def __init__(self, gui_main_wnd):
        QtCore.QThread.__init__(self)
        self.gui_main_wnd_ = gui_main_wnd

    def __del__(self):
        self.wait()

    def run(self):
        i = 0
        while True:
            i += 10
            time.sleep(0.3) # artificial time delay
            mouse_click_event = QtGui.QMouseEvent(QtCore.QEvent.MouseButtonPress, 
                                                  QtCore.QPoint(40, i+40),
                                                  QtCore.QPoint(40, i+40),
                                                  Qt.LeftButton, None)
            self.gui_main_wnd_.mousePressEvent(mouse_click_event)
            #self.emit( QtCore.SIGNAL('update(QString)'), "from work thread " + str(i) )

        self.terminate()

#=============================== END OF FILE =================================
