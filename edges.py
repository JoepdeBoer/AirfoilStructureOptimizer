import numpy as np
class Edge():
    def __init__(self, node1, node2, qbase, centroid):
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


    def momentccw(self, centroid):
        fx, fy = self.qbase * self.s * self.v # fx to the right fy up
        mz = -fx * self.ymid + fy * self.xmid # moment ccw +
        self.moment = float(mz)
