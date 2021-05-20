from PyQt4 import QtGui, QtCore
from GuiMainWindow import GuiMainWindow
from GuiVertex import GuiVertex
import time

#class GuiTrailExecutor(QtCore.QThread):
class GuiTrailExecutor(QtCore.QObject):
    def __init__(self, gui_main_wnd):
        #QtCore.QThread.__init__(self)
        super(GuiTrailExecutor, self).__init__() 
        self.gui_main_wnd_ = gui_main_wnd

    def __del__(self):
        self.wait()


    def run(self):
        i = 0
        while True:
            time.sleep(1) # artificial time delay
            i += 1
            # --- V1 -----------------------------------------------------------
            #i += 10
            #mouse_click_event = QtGui.QMouseEvent(QtCore.QEvent.MouseButtonPress,  
            #                                      QtCore.QPoint(40, i+40),
            #                                      QtCore.Qt.LeftButton, 
            #                                      QtCore.Qt.LeftButton,
            #                                      QtCore.Qt.NoModifier)
            #self.gui_main_wnd_.mousePressEvent(mouse_click_event)
            # -------------------------------------------------------------------
            gr_items = self.gui_main_wnd_.scene_.items()
            if len(gr_items) < 4:
                continue
            for ii in gr_items:
                if isinstance(ii,GuiVertex):
                    trgt_vrtx = ii

            trgt_pos = trgt_vrtx.pos()
            trgt_pos = QtCore.QPoint(trgt_pos.x() + 10*(0 if (i%4) == 0 or (i%4) == 3 else 1), 
                                     trgt_pos.y() + 10*(0 if (i%4) == 0 or (i%4) == 1 else 1))
            print trgt_pos
            mouse_mv_event = QtGui.QMouseEvent(QtCore.QEvent.MouseButtonPress, 
                                               trgt_pos,
                                               QtCore.Qt.LeftButton, 
                                               QtCore.Qt.LeftButton,
                                               QtCore.Qt.NoModifier)
            self.gui_main_wnd_.view_.viewport().mouseMoveEvent(mouse_mv_event)
            self.gui_main_wnd_.update()
            #trgt_vrtx.mouseMoveEvent(mouse_mv_event)
            #self.gui_main_wnd_.mouseMoveEvent(mouse_mv_event)

        self.terminate()

#=============================== END OF FILE =================================
