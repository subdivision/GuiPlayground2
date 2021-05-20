from PyQt4 import QtGui, QtCore
from GuiArrow import GuiArrow
from gc import enable

#-----------------------------------------------------------------------------
class GuiVertex(QtGui.QGraphicsEllipseItem):

    def __init__(self, vrtx_id, cntr_id, wnd, 
                 parent=None, 
                 pos = QtCore.QPointF(-5, -5),
                 size = 10.0, 
                 default_color = QtCore.Qt.blue,
                 hover_color = QtCore.Qt.cyan):
        self.vrtx_id_ = int(vrtx_id)
        self.cntr_id_ = int(cntr_id)
        self.wnd_ = wnd
        QtGui.QGraphicsPixmapItem.__init__(self, parent)
        QtGui.QGraphicsEllipseItem.__init__(self, 
                                            QtCore.QRectF(pos.x(), pos.y(), 
                                                          size, size), parent)

        self.reflectAction_ = QtGui.QAction('Reflect contour', None, 
                                            checkable=True, enabled=True)
        self.uniteReflAction_ = QtGui.QAction('Unite with reflection', None, 
                                              checkable=True, enabled=False)
        
        self.setFlag(QtGui.QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QtGui.QGraphicsItem.ItemIsSelectable, True)
        self.setAcceptHoverEvents(True)
        self.size_ = float(size)

        self.default_pen_ = QtGui.QPen(default_color,1)
        self.default_brush_ = QtGui.QBrush(default_color)
        self.hover_pen_ = QtGui.QPen(hover_color,1)
        self.hover_brush_ = QtGui.QBrush(hover_color)
        self.setDefaultSizeAndColor()
        self.arrow1_ = None
        self.arrow2_ = None
        self.mirr_vrtx_= None

    #-------------------------------------------------------------------------
    def setArrow1(self, arrow):
        self.arrow1_ = arrow

    #-------------------------------------------------------------------------
    def setArrow2(self, arrow):
        self.arrow2_ = arrow

    #-------------------------------------------------------------------------
    def setMirrorVertex(self, mirr_vrtx):
        self.mirr_vrtx_ = mirr_vrtx

    #-------------------------------------------------------------------------
    def contextMenuEvent(self, event):
        menu = QtGui.QMenu()
        snapAction = QtGui.QAction('Snap to grid', menu)
        snapAction.triggered.connect(self.snapToGrid)
        menu.addAction(snapAction)
        creaseAction = QtGui.QAction('Toggle crease', menu)
        creaseAction.triggered.connect(self.toggleCrease)
        menu.addAction(creaseAction)
        removeAction = QtGui.QAction('Remove vertex', menu)
        removeAction.triggered.connect(self.removeVertex)
        menu.addAction(removeAction)
        menu.addSeparator()
        removeCntrAction = QtGui.QAction('Remove this contour', menu)
        removeCntrAction.triggered.connect(self.removeContour)
        menu.addAction(removeCntrAction)
        splitEdgeAction = QtGui.QAction('Split next edge', menu)
        splitEdgeAction.triggered.connect(self.splitNextEdge)
        menu.addAction(splitEdgeAction)
        closeCntrAction = QtGui.QAction('Toggle close contour', menu)
        closeCntrAction.triggered.connect(self.toggleCloseContour)
        menu.addAction(closeCntrAction)
        naiveNormsAction = QtGui.QAction('Set naive normals', menu)
        naiveNormsAction.triggered.connect(self.setNaiveNormals)
        menu.addAction(naiveNormsAction)
        self.reflectAction_.setChecked(\
            self.wnd_.contours_[self.cntr_id_].isMirrored())
        self.uniteReflAction_.setChecked(\
            self.wnd_.contours_[self.cntr_id_].getMirrorConnected())
        menu.addAction(self.reflectAction_)
        if self.reflectAction_.isChecked():
            menu.addAction(self.uniteReflAction_)
        a = menu.exec_(event.screenPos())

        if a == self.reflectAction_:
            self.reflectContour(self.reflectAction_.isChecked())
        elif a == self.uniteReflAction_:
            self.wnd_.contours_[self.cntr_id_].setMirrorConnected(
                                             self.uniteReflAction_.isChecked())
            self.wnd_.update()

    #-------------------------------------------------------------------------
    def reflectContour(self, new_refl_state):
        if new_refl_state:
            self.uniteReflAction_.setEnabled(True)
            self.wnd_.contours_[self.cntr_id_].createMirrorContour()
        else:
            self.uniteReflAction_.setEnabled(False)
            self.wnd_.contours_[self.cntr_id_].removeMirrorContour()
        
    #-------------------------------------------------------------------------
    def snapToGrid(self):
        curr_x = int(self.pos().x())
        curr_x = (curr_x/10)*10 if 0 <= curr_x%10 <=4 else (curr_x/10 +1)*10
        curr_y = int(self.pos().y())
        curr_y = (curr_y/10)*10 if 0 <= curr_y%10 <=4 else (curr_y/10 +1)*10
        snapped_pos = QtCore.QPointF(curr_x, curr_y)
        self.setPos(snapped_pos)
        self.finishPositionUpdate()
        
    #-------------------------------------------------------------------------
    def toggleCrease(self):
        if not self.arrow2_:
            self.arrow2_ = GuiArrow(self, self.wnd_)
            self.arrow2_.setPos( self.arrow1_.pos())
            self.scene().addItem(self.arrow2_)
            if self.mirr_vrtx_:
                self.mirr_vrtx_.arrow2_ = GuiArrow(self.mirr_vrtx_, self.wnd_)
                self.mirr_vrtx_.arrow2_.setPos(self.mirr_vrtx_.arrow1_.pos())
                self.mirr_vrtx_.arrow2_, self.mirr_vrtx_.arrow1_ = \
                  self.mirr_vrtx_.arrow1_, self.mirr_vrtx_.arrow2_
        else:
            self.scene().removeItem(self.arrow2_)
            self.arrow2_ = None
            if self.mirr_vrtx_:
                self.mirr_vrtx_.arrow2_, self.mirr_vrtx_.arrow1_ = \
                  self.mirr_vrtx_.arrow1_, self.mirr_vrtx_.arrow2_
                self.mirr_vrtx_.arrow2_ = None
        self.wnd_.update()         
    
    #-------------------------------------------------------------------------
    def removeVertex(self):
        self.wnd_.removeVertex(self)
    
    #-------------------------------------------------------------------------
    def splitNextEdge(self):
        self.wnd_.splitNextEdge(self)

    #-------------------------------------------------------------------------
    def toggleCloseContour(self):
        self.wnd_.toggleCloseContour(self.cntr_id_)
    
    #-------------------------------------------------------------------------
    def removeContour(self):
        self.wnd_.removeContour(self.cntr_id_)
            
    #-------------------------------------------------------------------------
    def setNaiveNormals(self):
        self.wnd_.setNaiveNormals(self.cntr_id_)
        
    #-------------------------------------------------------------------------
    def mouseMoveEvent(self, event):
        QtGui.QGraphicsEllipseItem.mouseMoveEvent(self, event)
        target_pos = self.wnd_.view_.mapToScene(\
                        QtGui.QCursor().pos() - 
                        QtCore.QPoint(self.wnd_.geometry().x(), 
                                        self.wnd_.geometry().y()) -
                        self.wnd_.view_.pos())
        self.setPos(target_pos)
        self.finishPositionUpdate()
        self.wnd_.update()

    #-------------------------------------------------------------------------
    def finishPositionUpdate(self):
        curr_pos = self.pos()
        self.moveBy(-self.size_/2.0, -self.size_/2.0)
        if self.arrow1_:
            self.arrow1_.setPos( curr_pos )
            #self.arrow1_.moveBy(self.size_/2.0, self.size_/2.0)
        if self.arrow2_:
            self.arrow2_.setPos( curr_pos )
            #self.arrow2_.moveBy(self.size_/2.0, self.size_/2.0)
        str_prefix = 'ID: %d, Contour: %d, Pos: ' %(self.vrtx_id_,
                                                    self.cntr_id_)
        str_suffix = '[' + str('{0:.2f}'.format(curr_pos.x())) + ', ' \
                         + str('{0:.2f}'.format(curr_pos.y())) + ']'
        self.setToolTip( str_prefix + str_suffix ) 
        if self.mirr_vrtx_:
            self.mirr_vrtx_.setPos(-curr_pos.x() - self.size_/2.0, 
                                    curr_pos.y() - self.size_/2.0)

    #-------------------------------------------------------------------------
    def hoverEnterEvent(self, event):
        self.setHoverSizeAndColor()
        return QtGui.QGraphicsSceneHoverEvent.accept(event)

    #-------------------------------------------------------------------------
    def hoverLeaveEvent(self, event):
        self.setDefaultSizeAndColor()
        return QtGui.QGraphicsSceneHoverEvent.accept(event)

    #-------------------------------------------------------------------------
    def setDefaultSizeAndColor(self):
        self.setPen(self.default_pen_)
        self.setBrush(self.default_brush_)
        self.setRect(QtCore.QRectF(0, 0, self.size_, self.size_))

    #-------------------------------------------------------------------------
    def setHoverSizeAndColor(self):
        self.setPen(self.hover_pen_)
        self.setBrush(self.hover_brush_)
        self.setRect(QtCore.QRectF(0, 0, self.size_+1, self.size_+1))

    #-------------------------------------------------------------------------
    def saveIntoFile(self, fd):
        str_result = str(self.vrtx_id_) + ','
        curr_pos =  self.pos()
        curr_pos += QtCore.QPointF(self.size_/2.0, self.size_/2.0)
        str_result += str(curr_pos.x()) + ',' + str(curr_pos.y()) + ','
        str_result += str(self.arrow1_.dx_) + ','
        str_result += str(self.arrow1_.dy_) + ','
        str_result += str(self.arrow1_.gui_rotation_)
        if self.arrow2_:
            str_result += ',' + str(self.arrow2_.dx_) + ','
            str_result += str(self.arrow2_.dy_) + ','
            str_result += str(self.arrow1_.rotation())
        fd.write(str_result + '\n')

    #-------------------------------------------------------------------------
    @staticmethod
    def createVertexAndNormal(str_one_vrtx):
        pass
#=============================== END OF FILE =================================
