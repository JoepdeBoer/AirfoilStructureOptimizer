import numpy as np
from shapely import LinearRing, Polygon
from edges import Edge
class Cell2():
    def __init__(self, nodes: list, centroid: list, skin_thickness: float, G: float):
        self.structure_centroid = centroid
        self.skinthickness = skin_thickness
        self.G = G
        self.nodes = nodes
        self.edges = self.create_edges()
        self.polygon, self.ring = self.make_polygon()
        self.perimeter = self.ring.length
        self.area = self.polygon.area
        print('cell created with area', self.area)
        print('Total perimeter', self.perimeter)
        for edge in self.edges:
            print(edge.qbase)


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
        edges = []
        qb = 0
        prev = self.nodes[0]
        for node in self.nodes[1:]:
            edge = Edge(prev, node, qb, self.structure_centroid)
            edge.qbase += qb
            edges.append(edge)
            qb += node.dqb
            prev = node
        return edges




