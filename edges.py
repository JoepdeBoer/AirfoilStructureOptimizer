class section():
    def __init__(self, x1, y1, x2, y2, thickness):
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.s = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        self.Aeq = self.s * thickness # We assume 0 thickness skin after idealisation