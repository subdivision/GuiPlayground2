from PyQt4 import QtGui, QtCore
from GuiVertex import GuiVertex
from GuiArrow import GuiArrow

class GuiContour:
    def __init__(self, cntr_id):
        self.cntr_id_ =  cntr_id
        self.is_open_ = True
        self.real_cntr_ = []
        self.mirr_cntr_ = None
        self.mirr_cnctd_ = False

    def isMirrored(self):
        return self.mirr_cntr_ != None

    def createMirrorContour(self):
        self.mirr_cntr_ = []
        for real_vrtx in self.real_cntr_[::-1]:
            mirr_vrtx = self.createMirrorVertex(real_vrtx)
            self.mirr_cntr_.append(mirr_vrtx)

    def createMirrorVertex(self, real_vrtx):
        real_vrtx_x = real_vrtx.pos().x()
        if abs(abs(real_vrtx_x) - real_vrtx.size_/2.0) < 0.01:
            return real_vrtx

        mirr_vrtx = GuiVertex(-real_vrtx.vrtx_id_, 
                                real_vrtx.cntr_id_, 
                                real_vrtx.wnd_)
        mirr_arrw1 = GuiArrow(mirr_vrtx, real_vrtx.wnd_)
        mirr_arrw1.dx_ = -real_vrtx.arrow1_.dx_
        mirr_arrw1.dy_ = real_vrtx.arrow1_.dy_
        mirr_vrtx.setArrow1(mirr_arrw1)
        if real_vrtx.arrow2_:
            mirr_arrw2 = GuiArrow(mirr_vrtx, real_vrtx.wnd_)
            mirr_vrtx.setArrow2(mirr_arrw2)
            mirr_arrw2.dx_ = -real_vrtx.arrow2_.dx_
            mirr_arrw2.dy_ = real_vrtx.arrow2_.dy_
        mirr_pos = real_vrtx.pos()
        mirr_pos.setX(-mirr_pos.x() - mirr_vrtx.size_/2.0)
        mirr_pos.setY( mirr_pos.y() + mirr_vrtx.size_/2.0)
        mirr_vrtx.setPos(mirr_pos)
        mirr_vrtx.finishPositionUpdate()
        real_vrtx.mirr_vrtx_ = mirr_vrtx
        return mirr_vrtx

    def removeMirrorContour(self):
        for real_vrtx in self.real_cntr_[::-1]:
            real_vrtx.setMirrorVertex(None)
        self.mirr_cntr_ = None

    def getVertexByID(self, vrtx_id):
        vrtx = None
        cntr = self.real_cntr_ if vrtx_id > 0 else self.mirr_cntr_
        for v in cntr:
            if v.vrtx_id_ == vrtx_id:
                vrtx = v
                break
        return vrtx
    
    def append(self, vertex):
        self.real_cntr_.append(vertex)
        if self.isMirrored():
            mirr_vrtx = self.createMirrorVertex(vertex)
            self.mirr_cntr_.insert(0, mirr_vrtx)

    def insert(self, idx, vertex):
        self.real_cntr_.insert(idx, vertex)
        if self.isMirrored():
            mirr_vrtx = self.createMirrorVertex(vertex)
            if idx == 0:
                self.mirr_cntr_.append(mirr_vrtx)
            else:
                n = len(self.mirr_cntr_)
                self.mirr_cntr_.insert(n-idx, mirr_vrtx)
    
    def remove(self, vertex):
        mirr_vrtx = vertex.mirr_vrtx_
        if mirr_vrtx and self.mirr_cntr_:
            self.mirr_cntr_.remove(mirr_vrtx)
        self.real_cntr_.remove(vertex)

    def setMirrorConnected(self, is_conn):
        self.mirr_cnctd_ = is_conn

    def getMirrorConnected(self):
        return self.mirr_cnctd_

    def setOpen(self, b_is_open = True):
        self.is_open_ = b_is_open

    def getOpen(self):
        return self.is_open_

    def loadFromFile(self,fd):
        b_result = True
        try:
            self.cntr_id_ = int(fd.readline().strip())
            self.is_open_ = bool(fd.readline().strip())
            b_mirr = bool(fd.readline().strip())
            if b_mirr:
                self.mirr_cntr_ = []
            n_real_cntr = int(fd.readline().strip())
            for i in range(n_real_cntr):
                str_one_vrtx = fd.readline().strip()
                rv = GuiVertex.createVertexAndNormal(str_one_vrtx)
                self.real_cntr_.append(rv)
            n_mirr_cntr = int(fd.readline().strip())
            for i in range(n_mirr_cntr):
                str_one_vrtx = fd.readline().strip()
                mv = GuiVertex.createVertexAndNormal(str_one_vrtx)
                self.mirr_cntr_.append(mv)
            for mv in self.mirr_cntr_:
                for rv in self.real_cntr_[::-1]:
                    if mv.vrtx_id_ == -rv.vrtx_id_:
                        rv.mirr_vrtx_ = mv
        except Exception:
            b_result = False
        return b_result 
            
    def saveIntoFile(self, fd):
        fd.write(str(self.cntr_id_)+'\n')
        fd.write('1\n' if self.getOpen() else '0\n')
        fd.write('1\n' if self.isMirrored() else '0\n')
        fd.write(str(len(self.real_cntr_)) + '\n')
        for v in self.real_cntr_:
            v.saveIntoFile(fd)
        
#=============================== END OF FILE =================================
    