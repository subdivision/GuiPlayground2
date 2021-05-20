from PyQt4 import QtGui, QtCore
import math as m
from CSubd2D import get_angle
#-----------------------------------------------------------------------------
class GuiArrow(QtGui.QGraphicsPathItem):
    def __init__(self, anchor, wnd, parent=None, 
                 length = 30, 
                 default_color =QtCore.Qt.blue,
                 hover_color=QtCore.Qt.cyan):
        self.anchor_ = anchor
        self.wnd_ = wnd
        QtGui.QGraphicsPixmapItem.__init__(self, parent)
        self.setFlag(QtGui.QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QtGui.QGraphicsItem.ItemIsSelectable, True)
        self.setAcceptHoverEvents(True)
        self.length_ = length
        self.default_pen_ = QtGui.QPen(default_color, 2)
        self.hover_pen_ = QtGui.QPen(hover_color, 3)
        self.default_brush_ = QtGui.QBrush(default_color)
        self.hover_brush_ = QtGui.QBrush(hover_color)
        arrow_path = QtGui.QPainterPath()
        anchor_pos = QtCore.QPointF(0,0)
        arrow_path.moveTo(anchor_pos)
        arrow_path.lineTo(self.length_ * 0.7, 
                          anchor_pos.y())
        arrow_path.lineTo(self.length_ * 0.7, 
                          anchor_pos.y() + 3.0)
        arrow_path.lineTo(self.length_, 
                          anchor_pos.y())
        arrow_path.lineTo(self.length_ * 0.7, 
                          anchor_pos.y() - 3.0)
        arrow_path.lineTo(self.length_ * 0.7, 
                          anchor_pos.y())
        arrow_path.lineTo(anchor_pos)
        self.setPath(arrow_path)
        self.moveBy( self.anchor_.size_/2.0, 
                     self.anchor_.size_/2.0)
        self.setDefaultSizeAndColor()
        self.dx_ = 0.0
        self.dy_ = 1.0
        self.gui_rotation_ = 90.0
        self.setRotation(self.gui_rotation_)

    def mouseMoveEvent(self, event):
        scene_pos = event.scenePos()
        curr_pos = self.pos()
        x = scene_pos.x() - curr_pos.x()
        y = scene_pos.y() - curr_pos.y()
        leng = (x**2+y**2)**0.5
        self.setDir(x/leng, y/leng)
        self.wnd_.update()

    def setDir(self, x, y):
        self.dx_ = x
        self.dy_ = y
        new_alpha = get_angle(self.dx_, self.dy_)
        self.gui_rotation_ = new_alpha*180.0/m.pi
        self.setRotation(self.gui_rotation_)
        self.setToolTip( 'Vector: (%0.3f, %0.3f)' %(self.dx_, self.dy_))
        if self.anchor_.mirr_vrtx_:
            if self == self.anchor_.arrow1_:
                self.anchor_.mirr_vrtx_.arrow1_.dx_ = -self.dx_
                self.anchor_.mirr_vrtx_.arrow1_.dy_ = self.dy_
            else:
                self.anchor_.mirr_vrtx_.arrow2_.dx_ = -self.dx_
                self.anchor_.mirr_vrtx_.arrow2_.dy_ = self.dy_
        
    def hoverEnterEvent(self, event):
        self.setHoverSizeAndColor()
        return QtGui.QGraphicsSceneHoverEvent.accept(event)

    def hoverLeaveEvent(self, event):
        self.setDefaultSizeAndColor()
        return QtGui.QGraphicsSceneHoverEvent.accept(event)

    def setDefaultSizeAndColor(self):
        self.setPen(self.default_pen_)
        self.setBrush(self.default_brush_)

    def setHoverSizeAndColor(self):
        self.setPen(self.hover_pen_)
        self.setBrush(self.hover_brush_)


#=============================== END OF FILE =================================
