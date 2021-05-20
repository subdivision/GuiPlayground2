from PyQt4 import QtGui, QtCore
from GuiArrow import GuiArrow
from GuiVertex import GuiVertex
from GuiContour import GuiContour
import numpy as np
from functools import partial
import math as m
from CSubd2D import double_polygon, circle_avg, linear_avg, \
                  subd_4PT_one_step, subd_LR_one_step, \
                  init_normals, get_angle
from BSplineAvg import bspline_average_export_v3
#-----------------------------------------------------------------------------
class GuiIDGenerator:
    def __init__(self):
        self.next_id_to_create_ = 1
        
    def next_id(self):
        self.next_id_to_create_ += 1 
        return self.next_id_to_create_-1
    
    def curr_id(self):
        return self.next_id_to_create_-1
        
    
#-----------------------------------------------------------------------------
class GuiMainWindow(QtGui.QMainWindow):

    def __init__(self, parent=None):
        super(GuiMainWindow, self).__init__(parent)
        self.setWindowTitle('Circle Average Demo')
        self.setMouseTracking(True)
        self.scene_ = QtGui.QGraphicsScene(-300, -300, 800, 900)
        self.y_axis_ = QtGui.QGraphicsLineItem(0,800,0,-800)
        self.y_axis_.setPen(QtGui.QPen(QtCore.Qt.darkMagenta,1))
        self.scene_.addItem(self.y_axis_)
        #self.scene_.setBackgroundBrush(QtGui.QBrush(QtGui.QColor.fromRgb(245, 255, 250)))
        self.scene_.setBackgroundBrush(QtGui.QBrush(QtGui.QColor.fromRgb(255, 250, 250)))
         	
        self.view_ = QtGui.QGraphicsView()
        self.view_.setScene(self.scene_)
        self.view_.setGeometry(QtCore.QRect(0, 0, 800, 900))
        self.view_.setSceneRect(self.scene_.sceneRect())
        self.setCentralWidget(self.view_)
        self.showPointsAction_ = QtGui.QAction('Show points', None, 
                                            checkable=True, enabled=True)
        self.showPointsAction_.setChecked(True)
        self.showNormsAction_ = QtGui.QAction('Show normals', None, 
                                            checkable=True, enabled=True)
        self.showNormsAction_.setChecked(True)

        self.algo_label_ = QtGui.QLabel('Algorithm:')
        self.iter_label_ = QtGui.QLabel('Iterations:')
        self.mparam_label_ = QtGui.QLabel('m (*LR only):')
        self.normalsfix_label_ = QtGui.QLabel('Normals Fix:')

        self.iterations_ = QtGui.QComboBox()
        self.iterations_.addItems(['0','1','2','3','4', '5', '6'])
        self.iterations_.setCurrentIndex(1) #3
        self.connect(self.iterations_, 
                     QtCore.SIGNAL('currentIndexChanged(int)'),
                     self.update)  
        self.algo_ = QtGui.QComboBox()
        self.algo_.addItems(['LLR','MLR','L4Pt','M4Pt'])
        self.algo_.setCurrentIndex(1) #3
        self.connect(self.algo_, 
                     QtCore.SIGNAL('currentIndexChanged(int)'),
                     self.update)  
        self.mparam_ = QtGui.QComboBox()
        self.mparam_.addItems(['0','1','2','3'])
        self.mparam_.setCurrentIndex(0) #3
        self.connect(self.mparam_, 
                     QtCore.SIGNAL('currentIndexChanged(int)'),
                     self.update)  
        self.mparam_.setEnabled(False)
        self.nfix_ = QtGui.QCheckBox()
        self.nfix_.setChecked(False)
        self.nfix_.stateChanged.connect(self.update)

        #self.zoom_ = QtGui.QSlider(QtCore.Qt.Horizontal)
        #self.zoom_.setTickPosition(QtGui.QSlider.TicksBothSides)
        #self.zoom_.setRange(-3,3)
        #self.zoom_.setSliderPosition(0)
        #self.zoom_.setTickInterval(1)
        #self.zoom_.setFixedHeight(18)
        #self.connect(self.zoom_, 
        #             QtCore.SIGNAL('valueChanged(int)'),
        #             self.zoomView)  

        self.setStatusBar( QtGui.QStatusBar() )
        self.statusBar().showMessage(
                '[' + str('{0:.2f}'.format(self.pos().x())) + ', ' \
                    + str('{0:.2f}'.format(self.pos().y())) + ']')

        #self.statusBar().addPermanentWidget(self.zoom_)
        self.statusBar().addPermanentWidget(self.normalsfix_label_)
        self.statusBar().addPermanentWidget(self.nfix_)
        self.statusBar().addPermanentWidget(self.iter_label_)
        self.statusBar().addPermanentWidget(self.iterations_)
        self.statusBar().addPermanentWidget(self.algo_label_)
        self.statusBar().addPermanentWidget(self.algo_)
        self.statusBar().addPermanentWidget(self.mparam_label_)
        self.statusBar().addPermanentWidget(self.mparam_)
        self.scene_.installEventFilter(self)
        self.cntr_id_gen_ = GuiIDGenerator()
        self.vrtx_id_gen_ = GuiIDGenerator()
        self.contours_ = {}
        cntr_id = self.cntr_id_gen_.next_id()
        self.contours_[cntr_id] = GuiContour(cntr_id)
        self.pen_m4pt_ = QtGui.QPen(QtCore.Qt.red)
        self.pen_l4pt_ = QtGui.QPen(QtCore.Qt.green)
        self.pen_black_ = QtGui.QPen(QtCore.Qt.black)
        self.pen_black_.setWidth(2)
        self.pen_src_mesh_ = QtGui.QPen(QtCore.Qt.darkBlue)
        self.pen_src_mesh_.setWidth(1)
        self.pen_src_mesh_.setStyle(QtCore.Qt.DotLine)
        
    #-------------------------------------------------------------------------
    def zoomView(self, curr_z):
        gr_items = self.scene_.items()
        self.view_.resetMatrix()
        if 0.0 == curr_z:
            for ii in gr_items:
                if isinstance(ii,GuiArrow) or isinstance(ii,GuiVertex):
                    ii.setScale(1.0)
        else:
            curr_z = 1.0 + float(curr_z)/10.0
            self.view_.scale(curr_z, curr_z)
            for ii in gr_items:
                if isinstance(ii,GuiArrow) or isinstance(ii,GuiVertex):
                    ii.setScale(1.0/curr_z)


    #-------------------------------------------------------------------------
    def resizeEvent(self, event):
        self.view_.resizeEvent(event)
        nw = max(self.scene_.sceneRect().width(), event.size().width())
        nh = max(self.scene_.sceneRect().height(), event.size().height())
        x0 = -nw/2.0 
        y0 = -nh/2.0 
        self.scene_.setSceneRect(x0, y0, nw, nh)
        self.view_.setSceneRect(self.scene_.sceneRect())
        self.view_.centerOn(self.view_.mapFromScene(0,0))
        QtGui.QMainWindow.resizeEvent(self, event)

    #-------------------------------------------------------------------------
    def eventFilter(self, obj, event):
        if hasattr( event, "pos"):
            target_pos = self.view_.mapToScene(\
                           QtGui.QCursor().pos() - 
                           QtCore.QPoint(self.geometry().x(), 
                                         self.geometry().y()) -
                           self.view_.pos())
            #target_pos.setY(self.scene_.height() - target_pos.y())                
            self.statusBar().showMessage(
                '[' + str('{0:.2f}'.format(target_pos.x())) + ', ' \
                    + str('{0:.2f}'.format(target_pos.y())) + ']')
        return False

    #-------------------------------------------------------------------------
    def mousePressEvent(self, event):
        if event.button() != QtCore.Qt.LeftButton:
            return  
        cntr_id = self.cntr_id_gen_.curr_id()
        vrtx_id = self.vrtx_id_gen_.next_id()
        vrtx = GuiVertex(vrtx_id, cntr_id, self)
        arrw = GuiArrow(vrtx, self)
        #arrw.setDir(2.0**0.5/2.0, 2.0**0.5/2.0)
        vrtx.setArrow1(arrw)
        target_pos = self.view_.mapToScene(event.pos() - self.view_.pos())
        #target_pos.setY(self.scene_.height() - target_pos.y())                
        vrtx.setPos(target_pos)
        vrtx.finishPositionUpdate()
        self.contours_[cntr_id].append(vrtx)
        self.scene_.addItem(arrw)
        self.scene_.addItem(vrtx)
        arrw.setVisible(self.showNormsAction_.isChecked())
        vrtx.setVisible(self.showPointsAction_.isChecked())
        self.update()
    
    #-------------------------------------------------------------------------
    def removeVertex(self, vertex):
        self.contours_[vertex.cntr_id_].remove(vertex)
        gr_items = self.scene_.items()
        for ii in gr_items:
            if ii == vertex:
                if ii.arrow1_:
                    self.scene_.removeItem(ii.arrow1_)
                if ii.arrow2_:
                    self.scene_.removeItem(ii.arrow2_)
                self.scene_.removeItem(ii)
        self.update()

    #-------------------------------------------------------------------------
    def toggleCloseContour(self, cntr_id):
        self.contours_[cntr_id].setOpen(not self.contours_[cntr_id].getOpen())
        self.update()
        
    #-------------------------------------------------------------------------
    def removeContour(self, cntr_id):
        victim_cntr = self.contours_[cntr_id] 
        gr_items = self.scene_.items()
        for ii in gr_items:
            if ii in victim_cntr.real_cntr_:
                if ii.arrow1_:
                    self.scene_.removeItem(ii.arrow1_)
                if ii.arrow2_:
                    self.scene_.removeItem(ii.arrow2_)
                self.scene_.removeItem(ii)
        self.contours_.pop(cntr_id)
        if 0 == len(self.contours_.keys()):
            #We just removed the last contour
            cntr_id = self.cntr_id_gen_.next_id()
            self.contours_[cntr_id] = GuiContour(cntr_id)
        self.update()
    
    #-------------------------------------------------------------------------
    def setNaiveNormals(self, cntr_id):
        the_cntr = self.contours_[cntr_id] 
        b_open_cntr, pts, _, _, _, v2i, _ = \
                               self.collectContourWithMirror2Arrays(the_cntr)
        nrm = init_normals(pts, b_open_cntr)
        for i in range(len(v2i)):
            if v2i[i] < 0:
                continue
            curr_vrtx = the_cntr.getVertexByID(v2i[i])
            if i > 0 and v2i[i] == v2i[i-1]:
                curr_vrtx.arrow2_.setDir(nrm[i][0],nrm[i][1])
            else:
                curr_vrtx.arrow1_.setDir(nrm[i][0],nrm[i][1])
        self.update()
        
    #-------------------------------------------------------------------------
    def contextMenuEvent(self, event):    
        menu = QtGui.QMenu()
        menu.addAction(self.showPointsAction_)
        menu.addAction(self.showNormsAction_)
        addCntrAction = QtGui.QAction('Add contour', None)
        addCntrAction.triggered.connect(self.addContour)
        menu.addAction(addCntrAction)
        newDrawingAction = QtGui.QAction('New drawing', None)
        newDrawingAction.triggered.connect(self.newDrawing)
        menu.addAction(newDrawingAction)
        loadAction = QtGui.QAction('Load drawing', None)
        loadAction.triggered.connect(self.loadDrawing)
        menu.addAction(loadAction)
        saveAction = QtGui.QAction('Save drawing', None)
        saveAction.triggered.connect(self.saveDrawing)
        menu.addAction(saveAction)
        dumpAction = QtGui.QAction('Export PDF', None)
        dumpAction.triggered.connect(self.exportPDF)
        menu.addAction(dumpAction)
        quitAction = QtGui.QAction('Quit', None)
        quitAction.triggered.connect(self.close)
        menu.addAction(quitAction)
        a = menu.exec_(event.globalPos())
        if a == self.showPointsAction_:
            self.setShowPoints(self.showPointsAction_.isChecked())
        if a == self.showNormsAction_:
            self.setShowNorms(self.showNormsAction_.isChecked())


    #-------------------------------------------------------------------------
    def newDrawing(self):
        gr_items = self.scene_.items()
        for ii in gr_items:
            if ii == self.y_axis_:
                continue
            else:
                self.scene_.removeItem(ii)
        self.contours_ = {}
        cntr_id = self.cntr_id_gen_.next_id()
        self.contours_[cntr_id] = GuiContour(cntr_id)
        self.update()

    #-------------------------------------------------------------------------
    def addContour(self):
        cntr_id = self.cntr_id_gen_.next_id()
        self.contours_[cntr_id] = GuiContour(cntr_id)
    
    #-------------------------------------------------------------------------
    def loadDrawing(self):
        strFileName = QtGui.QFileDialog.getOpenFileName(self, "Load drawing",\
                                         "./", "Drawing Files (*.cntr *.txt)")        
        fd_drawing = open(strFileName, 'r')
        while True:
            curr_cntr = GuiContour(-1)
            b_status = curr_cntr.loadFromFile(fd_drawing)
            if not b_status:
                break
            self.contours_[curr_cntr.cntr_id_] = curr_cntr
        fd_drawing.close()
    
    #-------------------------------------------------------------------------
    def saveDrawing(self):
        strFileName = QtGui.QFileDialog.getSaveFileName(self, "Save drawing",\
                                         "./", "Drawing Files (*.cntr *.txt)")
        fd_drawing = open(strFileName, 'w')
        for curr_cntr in self.contours_.values():
                curr_cntr.saveIntoFile(fd_drawing)
        fd_drawing.close()

    #-------------------------------------------------------------------------
    def splitNextEdge(self, anchor):
        cntr_id = anchor.cntr_id_
        target_cntr = self.contours_[cntr_id].real_cntr_
        cntr_len = len(target_cntr)
        anchor_idx = target_cntr.index(anchor)
        if ( anchor_idx >= (cntr_len-1) and \
             self.contours_[cntr_id].getOpen() ) \
           or -1 == anchor_idx:
            return
        neighb = target_cntr[(anchor_idx+1)%cntr_len]
        anch_pos = anchor.pos()
        neig_pos = neighb.pos()
        s = anchor.size_/2.0
        p0 = np.array([anch_pos.x(), 
                       anch_pos.y()])
        if not anchor.arrow2_:
            n0 = np.array([anchor.arrow1_.dx_, 
                           anchor.arrow1_.dy_])
        else:
            n0 = np.array([anchor.arrow2_.dx_, 
                           anchor.arrow2_.dy_])

        p1 = np.array([neig_pos.x(), 
                       neig_pos.y()])
        n1 = np.array([neighb.arrow1_.dx_, 
                       neighb.arrow1_.dy_])
        vrtx_id = self.vrtx_id_gen_.next_id()
        vrtx = GuiVertex(vrtx_id, cntr_id, self)
        arrw = GuiArrow(vrtx, self)
        vrtx.setArrow1(arrw)
        ca_pt, ca_norm, _, _, _, _ = circle_avg(0.5, 0.5, True, 
                                                p0, p1, n0, n1)
        target_pos = QtCore.QPointF(ca_pt[0], ca_pt[1])
        vrtx.setPos(target_pos)
        arrw.setPos(target_pos + QtCore.QPointF(s,s))
        arrw.dx_ = ca_norm[0]
        arrw.dy_ = ca_norm[1]
        self.scene_.addItem(arrw)
        self.scene_.addItem(vrtx)
        self.contours_[cntr_id].insert((anchor_idx+1)%cntr_len, vrtx)
        self.update()

    #-------------------------------------------------------------------------
    def setShowPoints(self, b_show_points):
        gr_items = self.scene_.items()
        for ii in gr_items:
            if isinstance(ii,QtGui.QGraphicsEllipseItem) or\
               ii == self.y_axis_:
                ii.setVisible(b_show_points)

    #-------------------------------------------------------------------------
    def setShowNorms(self, b_show_norms):
        gr_items = self.scene_.items()
        for ii in gr_items:
            if isinstance(ii,QtGui.QGraphicsPathItem):
                ii.setVisible(b_show_norms)
        
    #-------------------------------------------------------------------------
    def paintEvent(self, event):
        QtGui.QMainWindow.paintEvent(self,event)
        gr_items = self.scene_.items()
        for ii in gr_items:
            if isinstance(ii,QtGui.QGraphicsLineItem) and ii != self.y_axis_:
                self.scene_.removeItem(ii)

        for curr_cntr in self.contours_.values():
            b_open_cntr, pts, _, pts2, _, _, _ = \
                               self.collectContourWithMirror2Arrays(curr_cntr)
            self.drawInScene(b_open_cntr, pts, self.pen_src_mesh_)
            self.drawInScene(b_open_cntr, pts2, self.pen_src_mesh_)

        for curr_cntr in self.contours_.values():
            b_open_cntr, pts, nrm, pts2, nrm2, _, _ = \
                               self.collectContourWithMirror2Arrays(curr_cntr)
            self.computeAndDraw(pts, nrm, b_open_cntr)
            self.computeAndDraw(pts2, nrm2, b_open_cntr)

    #-------------------------------------------------------------------------
    def collectContourWithMirror2Arrays(self, curr_cntr):
        b_open_cntr = curr_cntr.getOpen()
        pts = []
        nrm = []
        pts2 = []
        nrm2 = []
        v2i = []
        v2i2 = []

        self.collectGuiCntr2Array(curr_cntr.real_cntr_, pts, nrm, False, v2i)
        if curr_cntr.isMirrored():
            if b_open_cntr and curr_cntr.getMirrorConnected():
                b_open_cntr = False
                self.collectGuiCntr2Array(curr_cntr.mirr_cntr_, 
                                          pts, nrm, True, v2i)
            else:
                self.collectGuiCntr2Array(curr_cntr.mirr_cntr_, 
                                          pts2, nrm2, True, v2i2)
        return b_open_cntr, pts, nrm, pts2, nrm2, v2i, v2i2
    
    #-------------------------------------------------------------------------
    def collectGuiCntr2Array(self, contour, np_pts, np_norms, b_mirror, v2i):
        if len(contour) > 0:
            s = contour[0].size_/2.0
        for ii in contour:
            curr_x = ii.pos().x() + s
            curr_y = ii.pos().y() + s
            if np_pts and curr_x == np_pts[-1][0] and curr_y == np_pts[-1][1]:
                continue
            if np_pts and curr_x == np_pts[0][0] and curr_y == np_pts[0][1]:
                continue
            np_pts.append(np.array([curr_x, curr_y]))
            np_norms.append(np.array([ii.arrow1_.dx_, ii.arrow1_.dy_]))
            v2i.append(ii.vrtx_id_)
            if ii.arrow2_:
                np_pts.append(np.array([curr_x, curr_y]))
                np_norms.append(np.array([ii.arrow2_.dx_, ii.arrow2_.dy_]))
                v2i.append(ii.vrtx_id_)
                if b_mirror:
                    np_norms[-1], np_norms[-2] = np_norms[-2], np_norms[-1]

    #-------------------------------------------------------------------------
    def computeAndDraw(self, np_pts, np_norms, b_open_cntr):
        if len(np_pts)<2:
            return
        n_iterations = self.iterations_.currentIndex()
        str_algorithm = self.algo_.currentText()
        #fn_avg = circle_avg if str_algorithm[0] != 'L' else linear_avg
        fn_avg = bspline_average_export_v3 if str_algorithm[0] != 'L' else linear_avg
        if str_algorithm[1] == 'L':
            fn_algo = subd_LR_one_step  
            self.mparam_.setEnabled(True)
        else:
            fn_algo = subd_4PT_one_step
            self.mparam_.setEnabled(False)
        n_mparam = self.mparam_.currentIndex()
        
        for _ in range(n_iterations):
            np_pts, np_norms = fn_algo(np_pts, np_norms,b_open_cntr,
                                       fn_avg, n_mparam)
        if self.nfix_.isChecked():
            np_pts, np_norms = subd_LR_one_step(np_pts, np_norms,b_open_cntr,
                                                bspline_average_export, 0)
        self.drawInScene(b_open_cntr, np_pts, self.pen_black_)
        self.drawNormsInScene(np_pts, np_norms, self.pen_black_)

        #---------------------------------------------------------------------
        ##subd_func = subd_LR_one_step
        #subd_func = subd_4PT_one_step
        #cir_pts, cir_norms = np_pts[:], np_norms[:]
        #lin_pts, lin_norms = np_pts[:], np_norms[:]
        #for _ in range(n_iterations):
        #    if 0 == len(cir_pts) or 0 == len(lin_pts):
        #        return
        #    cir_pts, cir_norms = subd_func(cir_pts, cir_norms,  
        #                                           b_open_cntr, circle_avg)
        #    lin_pts, lin_norms = subd_func(lin_pts, lin_norms, 
        #                                           b_open_cntr, linear_avg)
        #self.drawInScene(b_open_cntr, cir_pts, self.pen_m4pt_)
        #self.drawInScene(b_open_cntr, lin_pts, self.pen_l4pt_)
        #---------------------------------------------------------------------

    #-------------------------------------------------------------------------
    def drawInScene(self, b_open_cntr, np_pts, pen ):
        nn = len(np_pts)-1 if b_open_cntr else len(np_pts)
        n = len(np_pts)
        for i in range(nn):
                segm = QtGui.QGraphicsLineItem()
                segm.setLine(np_pts[i][0], np_pts[i][1], 
                             np_pts[(i+1)%n][0], np_pts[(i+1)%n][1])
                segm.setPen(pen)
                self.scene_.addItem(segm)

    #-------------------------------------------------------------------------
    def drawNormsInScene(self, np_pts, np_norms, pen ):
        n = len(np_pts)
        for i in range(n):
                segm = QtGui.QGraphicsLineItem()
                segm.setLine(np_pts[i][0], np_pts[i][1], 
                             np_pts[i][0] + 4*np_norms[i][0], 
                             np_pts[i][1] + 4*np_norms[i][1])
                segm.setPen(pen)
                self.scene_.addItem(segm)

    #-------------------------------------------------------------------------
    def setSymmetry(self, new_symm_state):
        if new_symm_state == QtCore.Qt.Checked:
            self.ch_unite_mirr_.setEnabled(True)
            for curr_cntr in self.contours_.values():
                curr_cntr.createMirrorContour()
        else:
            self.ch_unite_mirr_.setEnabled(False)
            for curr_cntr in self.contours_.values():
                curr_cntr.removeMirrorContour()

    #-------------------------------------------------------------------------
    def exportPDF(self):
        file_name = 'C:\\TAU\\CircleAvg\\test.pdf'
        printer = QtGui.QPrinter(QtGui.QPrinter.HighResolution)
        printer.setPageSize(QtGui.QPrinter.A4)
        printer.setOrientation(QtGui.QPrinter.Portrait)
        printer.setOutputFormat(QtGui.QPrinter.PdfFormat)
        printer.setOutputFileName(file_name)

        p = QtGui.QPainter(printer)
        self.scene_.render(p)
        p.end()

#=============================== END OF FILE =================================
