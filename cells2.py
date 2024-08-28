import matplotlib.pyplot as plt
import numpy as np
from shapely import LinearRing, Polygon, plotting
from edges import Edge
class Cell2():
    def __init__(self, nodes: list, centroid: list, skin_thickness: float, G: float, number: int = 0):
        self.structure_centroid = centroid
        self.skinthickness = skin_thickness
        self.G = G
        self.nodes = nodes
        self.edges = self.create_edges()
        self.polygon, self.ring = self.make_polygon()
        self.perimeter = self.ring.length
        self.area = self.polygon.area
        self.number = number


        # self.spar_perimeter = [i[2] - i[1] for i in self.spars]
        # self.skinperimeter = self.perimeter - np.sum(self.spar_perimeter)
        # print('cell created with area', self.area)
        # print('Total perimeter', self.perimeter)
        # print('Spar perimeter', self.spar_perimeter)
        # print('Skin perimeter', self.skinperimeter)


    def make_polygon(self):
        x = []
        y = []
        for node in self.nodes:
            x.append(node.x)
            y.append(node.y)
        coordinates = np.stack((x, y), axis=1)
        return Polygon(coordinates), LinearRing(coordinates)


    def create_edges(self):
        # getting index of first spar for cell 1 at the bottom and cell 2 at the top
        for index, node in enumerate(self.nodes):
            if index != 0 and node.neighbors[2]:
                break

        # creating edges
        edges = []
        qblist = [0]
        qb = 0
        prev = self.nodes[index - 1]
        for node in self.nodes[index:]:
            edge = Edge(prev, node, qb, self.structure_centroid)
            edge.qbase += qb
            edge.is_spar = 1
            edges.append(edge)
            # print(f'{qb}  n1 - n2 = {prev.number}  {node.number} delta = {node.dqb}')
            qb += node.dqb
            qblist.append(qb)
            prev = node
        # print(f'{qb}  n1 - n2 = {prev.number}  {self.nodes[0].number} delta = {node.dqb}')
        for node in self.nodes[:index]:
            edge = Edge(prev, node, qb, self.structure_centroid)
            edge.qbase += qb
            edges.append(edge)
            # print(f'{qb}  n1 - n2 = {prev.number}  {node.number} delta = {node.dqb}')
            qb += node.dqb
            qblist.append(qb)
            prev = node





        # lastedge = Edge(prev, self.nodes[0], qb, self.structure_centroid)
        # edges.append(lastedge)
        # plt.scatter(range(len(qblist)), qblist)
        # plt.show()
        return edges




