import numpy as np
class Node():
    def __init__(self, x, y, number):
        self.number = number
        self.x = x # x location from centroid
        self.y = y # y location from centroid
        self.sigmaz = 0.
        self.A = 0.  # equavalent Area of the node due to skin
        self.dqb = 0. # jump in shear flow over this node due to shear loading
        self.neighbors = [None, None, None] # prev, next, connected spar ccw
        # self.starting_node = False # first node in cell
        # self.ending_node = False # last node in cell
    def compute_dqb(self, Ixx, Iyy, Ixy, Vx, Vy):
        dq1 = - (Vy * Iyy - Vx * Ixy) / (Ixx * Iyy - Ixy ** 2) * self.A * self.y
        dq2 = - (Vx * Ixx - Vy * Ixy) / (Ixx * Iyy - Ixy ** 2) * self.A * self.x
        self.dqb = dq1 + dq2

    def compute_A(self, thickness):
        for i in self.neighbors:
            if i:
                ds = np.sqrt((i.x - self.x) ** 2 + (i.y - self.y) ** 2)
                self.A += thickness * ds/6 * (2 + (i.sigmaz/self.sigmaz))

