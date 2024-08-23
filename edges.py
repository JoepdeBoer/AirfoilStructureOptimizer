import numpy as np
class Edge():
    def __init__(self, node1, node2, qbase, centroid):
        self.node1 = node1
        self.node2 = node2
        self.x1 = node1.x
        self.y1 = node1.y
        self.x2 = node2.x
        self.y2 = node2.y
        self.s = np.sqrt((self.x2 - self.x1) ** 2 + (self.y2 - self.y1) ** 2)
        self.v = np.array([self.x1 - self.x2, self.y1 - self.y2])/self.s
        self.xmid = (self.x1 + self.x2) / 2
        self.ymid = (self.y1 + self.y2) / 2
        self.qbase = qbase
        self.moment = self.momentccw(centroid)
        self.qs0 = 0
        self.qtotal = None
        self.stress = None
        self.is_spar = 0 # 0 is skin, 1 is spar1, 2 is spar2


    def momentccw(self, centroid):
        fx, fy = self.qbase * self.s * self.v # fx to the right fy up
        mz = -fx * self.ymid + fy * self.xmid # moment ccw +
        return float(mz)
