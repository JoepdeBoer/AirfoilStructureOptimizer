import numpy as np
import importAirfoil.importAirfoil as importAirfoil
import matplotlib.pyplot as plt
from shapely import LinearRing, Polygon, LineString, plotting
from nodes import Node
from cells2 import Cell2


class Airfoil:
    sparresolution = .03# [m]resolution
    """ creates idealises thin walled structure"
         from lednicer formatted airfoil cooridinates"""

    def __init__(self, filename, thickness, spars, chord, Load,):
        """
        :param filename: name of file with airfoil coordinates
        :param thickness: thickness of skin in [mm]
        :param spars: list of relative spar coordinates x/chord [-]
        :param chord: chord length in [m]
        :param Load: [Vx, Vy, Tz, Mx] shear postive up and right Torque positive ccw [N, N, Nm] forces applied at (0,0)
        """

        # Due to scipy not adhering to bounds aproximating the derivatives
        if spars[0] < 0.0001:
            spars[0][0] = 0.0001
        if spars[0] > spars[1]:
            print("spars[0] must be less than spars[1]")
            if spars[0] + 0.01 >= 1:
                print("spars[0] must be less than 1")
                spars[0] = 0.98
                spars[1] = 0.99
            else:
                spars[1] = spars[0] + 0.01



        self.G_skin = 27 * 10 ** 9  # [Pa] shear modulus
        self.Ixx = 0.
        self.Iyy = 0.
        self.Ixy = 0.
        self.file = filename
        self.Load = Load
        self.thickness = thickness
        self.spar_thickness = thickness  # TODO make spar thickness independent of skin-thickness
        self.chord = chord
        self.x, self.y = self.loadpoints() * chord
        self.spars = self.createspars(spars)
        self.indx_lastspar = int(self.spars[-1][-1]) + 1
        self.structural_skin_x = self.x[-self.indx_lastspar:self.indx_lastspar] # assuming last spar is where aileron begins
        self.structural_skin_y = self.y[-self.indx_lastspar:self.indx_lastspar]
        self.perimeter = LineString(np.column_stack((self.structural_skin_x, self.structural_skin_y)))
        self.centroid = self.calculate_centriod()
        self.calculate_I()
        self.skin_nodes = []
        self.spar_nodes = []
        self.create_nodes()
        self.idealise()
        self.cells = self.create_cells2()



    def loadpoints(self):
        """ loads airfoil coordinates from file"""
        points = np.array(importAirfoil(self.file))
        return points

    def createspars(self, spars):
        """"
        finding the nodes closesst to the location of the spars
        and linear interpolate to get the ymin and ymax of the spars

        """
        spargeom = np.zeros((len(spars), 4))  # [x, ymin, ymax, indx1]
        for i, x in enumerate(spars):
            xscaled = x * self.chord
            indx1 = np.searchsorted(self.x[len(self.x) // 2:], xscaled) + len(self.x) // 2
            indx2 = indx1 - 1
            indx3 = -indx1
            indx4 = indx3 - 1
            x1, y1 = self.x[indx1], self.y[indx1]
            x2, y2 = self.x[indx2], self.y[indx2]
            x3, y3 = self.x[indx3], self.y[indx3]
            x4, y4 = self.x[indx4], self.y[indx4]
            a1 = (y2 - y1) / (x2 - x1)
            a2 = (y4 - y3) / (x4 - x3)
            b1 = y1 - a1 * x1
            b2 = y3 - a2 * x3

            ymin = a1 * xscaled + b1 # linear interpolation to find lower location of spar
            ymax = a2 * xscaled + b2

            # shifting the points on the airfoil to match spar position
            self.x[indx1] = xscaled
            self.x[indx4] = xscaled
            self.y[indx1] = ymin
            self.y[indx4] = ymax
            spargeom[i] = [xscaled, ymin, ymax, indx1]
        return spargeom

    # def createcells(self):
    #     cells = {}
    #     for i, x in enumerate(self.spars):
    #         indx = int(x[3])
    #         cellcount = 1
    #         if i == 0:
    #             #first cell
    #             dsection_y = self.y[-indx - 1:indx + 1]
    #             dsection_x = self.x[-indx - 1:indx + 1]
    #             skin = np.column_stack((dsection_x, dsection_y))
    #             cells[f'cell_{cellcount}'] = Cell([skin], [x], [self.thickness], [self.thickness]) # TODO add diffrent thicknesses
    #             cellcount += 1
    #             plt.plot(dsection_x, dsection_y, color='pink')
    #             plt.plot([x[0], x[0]], [x[1], x[2]], color='blue')
    #             plt.show()
    #
    #         if i < len(self.spars) - 1:
    #             #middle cells
    #             spar2 = self.spars[i + 1]
    #             indx2 = int(spar2[3])
    #             bottom = np.column_stack((self.x[indx:indx2 + 1], self.y[indx:indx2 + 1]))
    #             top = np.column_stack((self.x[-indx2 - 1:-indx], self.y[-indx2 - 1:-indx]))
    #             cells[f'cell_{cellcount}'] = Cell([top, bottom], [x, spar2], [self.thickness, self.thickness], [self.thickness, self.thickness])
    #             cellcount += 1
    #             plt.plot(top[:, 0], top[:, 1], color='pink')
    #             plt.plot([x[0], x[0]], [x[1], x[2]], color='blue')
    #             plt.plot(bottom[:, 0], bottom[:, 1], color='pink')
    #             plt.plot([spar2[0], spar2[0]], [spar2[1], spar2[2]], color='red')
    #             plt.show()
    #     # for i in cells.values():
    #     #     i.make_sections(self.sparresolution)
    #
    #     return cells

    def create_cells2(self):
        cells = dict()
        indx1 = int(self.spars[0][3]) - (len(self.x) - len(self.skin_nodes))//2
        indx2 = int(self.spars[-1][3]) - (len(self.x) - len(self.skin_nodes))//2
        # first cell
        skin_nodes = self.skin_nodes[-indx1 -1: indx1 +1]
        spar_nodes = self.spar_nodes[0]
        nodes =  skin_nodes + spar_nodes
        cells[f'cell_1'] = Cell2(nodes, self.centroid, self.thickness, self.G_skin, number=1)

        # middle cell
        top_skin_nodes = self.skin_nodes[0: -indx1]
        bottom_skin_nodes = self.skin_nodes[indx1:]
        nodes = top_skin_nodes + self.spar_nodes[0][::-1] + bottom_skin_nodes + self.spar_nodes[1] # fliped direction of spar for ccw direction
        cells[f'cell_2'] = Cell2(nodes, self.centroid, self.thickness, self.G_skin, number=2)


        # x = []
        # y = []
        # c = []
        # for node in cells['cell_2'].nodes:
        #     print(node.number)
        #     if node in top_skin_nodes:
        #         c.append('pink')
        #     elif node in bottom_skin_nodes:
        #         c.append('blue')
        #     else:
        #         print(f'{node.x}, {node.y}')
        #         c.append('red')
        #     x.append(node.x)
        #     y.append(node.y)
        # plt.scatter(x, y, color=c)
        # plt.show()
        return cells




    def calc_qt(self):
        """
        caclculate the shear flow due to Torsion (constant per cell)
        shear flow in spar is addition of neighbouring cells
        """
        A = np.zeros((len(self.cells) + 1, len(self.cells) + 1))
        x = np.zeros((len(self.cells) + 1, 1))
        x[-1] = self.Load[2]  # Torque
        A[-1:-1] = 0
        spar1perimeter = self.spars[0][2] - self.spars[0][1]

        A[0, 0] = 1 / 2 / self.cells['cell_1'].area / self.thickness / self.G_skin * self.cells['cell_1'].perimeter
        A[0, 1] = -1 / 2 / self.cells['cell_1'].area / self.thickness / self.G_skin * spar1perimeter
        A[0, 2] = -1
        A[1, 0] = -1 / 2 / self.cells['cell_2'].area / self.thickness / self.G_skin * spar1perimeter
        A[1, 1] = 1 / 2 / self.cells['cell_2'].area / self.thickness / self.G_skin * self.cells['cell_2'].perimeter
        A[1, 2] = -1
        A[2, 0] = 2 * self.cells['cell_1'].area
        A[2, 1] = 2 * self.cells['cell_2'].area
        A[2, 2] = 0

        qt = np.linalg.solve(A, x)  # last element is dtheta/dz
        return qt

    def create_nodes(self):
        n = 0
        sparconnected = []
        x_c = self.structural_skin_x - self.centroid[0]  # coordinates w.r.t. centroid
        y_c = self.structural_skin_y - self.centroid[1]
        for i, x in enumerate(x_c):
            node = Node(x, y_c[i], n)
            if self.skin_nodes:
                node.neighbors[0] = self.skin_nodes[i-1]
                self.skin_nodes[i-1].neighbors[1] = node
            self.skin_nodes.append(node)
            if self.structural_skin_x[i] in self.spars[:,0]:
                sparconnected.append(node)
            n += 1

        for i, spar in enumerate(self.spars):
            sparlist = []
            xspar = spar[0] - self.centroid[0]
            ytop = spar[2]  - self.centroid[1] - self.sparresolution
            ybot = spar[1] - self.centroid[1] + self.sparresolution
            number = int((ytop-ybot)//self.sparresolution + 1)
            if number < 2:
                number = 2
            ys = np.linspace(ybot, ytop, number, endpoint=True)
            for j, y in enumerate(ys):
                node = Node(xspar, y, n)
                if sparlist and y != ybot:
                    node.neighbors[0] = sparlist[i-1]
                    sparlist[i-1].neighbors[1] = node
                sparlist.append(node)
                n += 1

            self.spar_nodes.append(sparlist)

        self.n = n
        sparconnected[0].neighbors[2] = self.spar_nodes[1][-1]
        sparconnected[1].neighbors[2] = self.spar_nodes[0][-1]
        sparconnected[2].neighbors[2] = self.spar_nodes[0][0]
        sparconnected[3].neighbors[2] = self.spar_nodes[1][0]
        self.spar_nodes[1][-1].neighbors[1] = sparconnected[0]
        self.spar_nodes[1][0].neighbors[0] = sparconnected[3]
        self.spar_nodes[0][-1].neighbors[1] = sparconnected[1]
        self.spar_nodes[0][0].neighbors[0] = sparconnected[2]

        # plt.scatter(self.x- self.centroid[0], self.y - self.centroid[1], c='black')
        # plt.scatter(sparconnected[0].x, sparconnected[0].y, c='purple')
        # plt.scatter(sparconnected[1].x, sparconnected[1].y, c='pink')
        # plt.scatter(sparconnected[2].x, sparconnected[2].y, c='red')
        # plt.scatter(sparconnected[3].x, sparconnected[3].y, c='orange')
        # plt.scatter(self.spar_nodes[1][-1].x, self.spar_nodes[1][-1].y, c='purple')
        # plt.scatter(self.spar_nodes[0][-1].x, self.spar_nodes[0][-1].y, c='pink')
        # plt.scatter(self.spar_nodes[0][0].x, self.spar_nodes[0][0].y, c='red')
        # plt.scatter(self.spar_nodes[1][0].x, self.spar_nodes[1][0].y, c='orange')
        # plt.axis('equal')
        # plt.show()

    def idealise(self):
        '''
        Compute the skin contribution to node areas based on the normal stress
        And calculate the shear flow jump over the nodes
        '''

        for node in self.skin_nodes:
            # node.sigmaz = self.Load[3] * self.Iyy * node.y - self.Load[3] * self.Ixy * node.x /(self.Ixx * self.Iyy - self.Ixy ** 2)
            node.sigmaz = self.Load[3] * node.y / self.Ixx
        for sparlist in  self.spar_nodes:
            for node in sparlist:
                # node.sigmaz = self.Load[3] * self.Iyy * node.y - self.Load[3] * self.Ixy * node.x /(self.Ixx * self.Iyy - self.Ixy ** 2)
                node.sigmaz = self.Load[3] * node.y / self.Ixx
        for node in self.skin_nodes:
            node.compute_A(self.thickness)
            node.compute_dqb(self.Ixx, self.Iyy, self.Ixy, self.Load[0], self.Load[1])
            # if node.A < 0:
            #     print("skin nodes" , node.number, node.A, node.neighbors)
        for sparlist in self.spar_nodes:
            for node in sparlist:
                node.compute_A(self.thickness) # TODO Could vary spar thickness
                node.compute_dqb(self.Ixx, self.Iyy, self.Ixy, self.Load[0], self.Load[1])
                # if node.A < 0:
                #     print("spar nodes", node.number, node.A, node.dqb, node.x+self.centroid[0], node.y)

                # print('spar nodes', node.A, node.dqb)



    def calculate_I(self):
        """calculate Ixx, Iyy, Ixy, of the section
        """
        xs = self.structural_skin_x - self.centroid[0]
        ys = self.structural_skin_y - self.centroid[1]
        for i, _ in enumerate(xs):
            x1 = xs[i]
            y1 = ys[i]
            if i < len(xs) - 1:
                x2 = xs[i+1]
                y2 = ys[i+1]
                yc = (ys[i] + ys[i+1]) / 2
                xc = (xs[i] + xs[i+1]) / 2
            else:
                x2 = xs[0]
                y2 = ys[0]
                yc = (ys[i] + ys[0]) / 2
                xc = (xs[i] + xs[0]) / 2
            ds = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
            dA = ds * self.thickness
            beta = np.arctan2(y1 - y2, x1 - x2)
            self.Ixx += ds**3 * self.thickness * np.sin(beta)**2 / 12 + dA * yc**2
            self.Iyy += ds**3 * self.thickness * np.cos(beta)**2 / 12 + dA * xc**2
            self.Ixy += ds**3 * self.thickness * np.sin(beta) * np.cos(beta) / 12 + dA * xc * yc
        for i in self.spars:
            y1 = i[1] - self.centroid[1]
            y2 = i[2] - self.centroid[1]
            x = i[0] - self.centroid[0]
            yc = (y1 + y2) / 2
            ds = y2 - y1
            dA = ds * self.spar_thickness
            self.Ixx += self.spar_thickness * ds**3 /12 + dA * yc**2
            self.Iyy += self.spar_thickness * ds**3 /12 + dA * x**2
            self.Ixy += self.spar_thickness * ds**3 /12 + dA * x * yc

        # As code originally used for arbitrary shapes Ixy computed however
        # it is not needed for the optimisation problem
        self.Ixy = 0 # !WARNING ONLY FOR SYMETRICAL STRUCTURES!

    def calculate_centriod(self):
        """
        calculate the centroid of the structural section bounded by the last spar
        :return:
        """
        xy = self.perimeter.centroid.coords.xy
        skincentroid = np.array([xy[0][0], xy[1][0]])
        skinarea = self.perimeter.length * self.thickness
        spar_centroids = np.zeros((len(self.spars) - 1, 3))
        for i, s in enumerate(spar_centroids): # last spar is already calculated in skin centroid
            spar_lenght = s[2] - s[1]
            spar_y = s[1] + spar_lenght / 2
            spar_area = spar_lenght * self.spar_thickness
            spar_centroids[i] = [s[0], spar_y, spar_area]
        centroid = (skincentroid * skinarea + np.sum(spar_centroids[:, 0:2] * spar_centroids[:, 2])) / (
                skinarea + np.sum(spar_centroids[:, 2]))

        return centroid

    def solve_shear(self):
        # based on dz/dthetha is same for each cell (no warping)
        # and Moment equavalance of shear forces with internal shearflow
        Amat = np.zeros((3, 3))
        Amat[:-1, -1] = -1
        A1 = self.cells['cell_1'].area
        A2 = self.cells['cell_2'].area

        I1 = 0
        I2 = 0
        Mx1 = 0
        Mx2 = 0
        # Cell1 circle integral of qb / tG ds
        for edge in self.cells['cell_1'].edges:
            I1 += edge.qbase * edge.s / self.thickness / self.G_skin
        # Cell2 circle integral of qb / tG ds
        for edge in self.cells['cell_2'].edges:
            I2 += edge.qbase * edge.s / self.thickness / self.G_skin
        I1 /= 2 * A1
        I2 /= 2 * A2


        # lengths of sections
        Lspar1 = self.spars[0][2] - self.spars[0][1]
        Lspar2 = self.spars[1][2] - self.spars[1][1]
        Lskin1 = self.cells['cell_1'].perimeter - Lspar1
        Lskin2 = self.cells['cell_2'].perimeter - Lspar2 - Lspar1

        Amat[0, :-1] = np.array([Lskin1 + Lspar1, -Lspar1])/A1/2/self.thickness/self.G_skin
        Amat[1, :-1] = np.array([-Lspar1, Lskin2 + Lspar2 + Lspar1 ])/A2/2/self.thickness/self.G_skin
        Amat[2, :-1] = 2*A1, 2*A2

        # moments due to base shearflow
        for edge in self.cells['cell_1'].edges:
            Mx1 += edge.moment
        skip = False
        firstspar = True
        for i, edge in enumerate(self.cells['cell_2'].edges):
            if edge.node2.neighbors[2] and skip==False and firstspar:
                Mx2 += edge.moment
                skip = True
            elif edge.node2.neighbors[2] and firstspar:
                skip = False
                firstspar = False
                continue
            if skip:
                continue
            Mx2 += edge.moment

        shear_induced_moment = self.Load[0] * self.centroid[1] - self.Load[1] * self.centroid[0] # ccw positive
        Mbase = Mx1 + Mx2
        b = np.array([-I1, -I2, shear_induced_moment - Mbase])

        # solve for qs dthetha/dz
        qs = np.linalg.solve(Amat, b)
        qs1 = qs[0]
        qs2 = qs[1]
        dthetha = qs[2]
        spar1 = False
        for edge in self.cells['cell_1'].edges:
            edge.qs0 += qs1
            if edge.is_spar == 1:
                edge.qs0 -= qs2
                edge.qtotal = edge.qs0 + edge.qbase
                # print(
                #     f'{edge.node1.number} - {edge.node2.number} qb = {edge.qbase} qs = {edge.qs0} qtotal = {edge.qtotal}')
            else:
                edge.qtotal = edge.qs0 + edge.qbase
        for edge in self.cells['cell_2'].edges:
            edge.qs0 += qs2
            if edge.is_spar == 1:
                edge.qs0 -= qs1
                edge.qtotal = edge.qs0 + edge.qbase
                # print(
                #     f'{edge.node1.number} - {edge.node2.number} qb = {edge.qbase} qs = {edge.qs0} qtotal = {edge.qtotal}')
            else:
                edge.qtotal = edge.qs0 + edge.qbase
    def find_max_shear_flow(self):
        self.solve_shear()
        qt1, qt2, _ = self.calc_qt()

        max_shear_flow = 0
        limiting_edge = None
        limiting_cell = None
        for edge in self.cells['cell_1'].edges:
            edge.qtotal += qt1
            if edge.is_spar == 1:
                edge.qtotal -= qt2
            if abs(edge.qtotal) > abs(max_shear_flow):
                max_shear_flow = edge.qtotal
                limiting_edge = edge
                limiting_cell = 1
        for edge in self.cells['cell_2'].edges:
            edge.qtotal += qt2
            if edge.is_spar == 1:
                edge.qtotal -= qt1
            if abs(edge.qtotal) > abs(max_shear_flow):
                max_shear_flow = edge.qtotal
                limiting_edge = edge
                limiting_cell = 2

        # print(f'max_shear_flow = {max_shear_flow}')
        # print(f'limiting_edge = {limiting_edge.node1.number} - {limiting_edge.node2.number} cell = {limiting_cell}')
        return max_shear_flow

    def structural_area(self):
        Snew = self.perimeter.length + self.spars[0][2] - self.spars[0][1] + self.spars[1][2] - self.spars[1][1]
        A = Snew * self.thickness
        # print(f'spar2 pos = {self.spars[1][0]}\n'
        #       f'Structural area = {A}\n'
        #       f'skin lenght = {self.perimeter.length}\n'
        #       f' spar1 lenght = {self.spars[0][2] - self.spars[0][1]}\n'
        #       f' spar2 lenght = {self.spars[1][2] - self.spars[1][1]}\n')
        return A





if __name__ == "__main__":
    # airfoil = Airfoil("NACA012.txt", 0.01, [0.4, 0.8], 2.5, [0, 75*10**3, 200*10**3, 1])
    airfoil = Airfoil("NACA0012-1000.txt", 0.005, [0.5, 0.9], 2.5, [0, 75 * 10 ** 3, 200 * 10 ** 3, 1])
    print(f'Area = {airfoil.structural_area()} \n'
          f'Stress = {airfoil.find_max_shear_flow()/airfoil.thickness/10**6} MPa\n'
          f'Shearstrenght aluminium = 207 MPa')
    plt.show()
