import numpy as np
import importAirfoil.importAirfoil as importAirfoil
import matplotlib.pyplot as plt
from shapely import LinearRing, Polygon, LineString


class Airfoil:
    sparresolution = 1e-3  # [m]resolution
    """ creates idealises thin walled structure"
         from lednicer formatted airfoil cooridinates"""

    def __init__(self, filename, thickness, spars, chord, Load, n=None):
        """
        :param filename: name of file with airfoil coordinates
        :param thickness: thickness of skin in [mm]
        :param spars: list of relative spar coordinates x/chord [-]
        :param chord: chord length in [m]
        :param Load: [Vx, Vy, Tz] shear postive up and right Torque positive ccw [N, N, Nm]
        :param n: number of nodes for per surface (top and bottom) total number is thus 2n
        """
        self.file = filename
        self.thickness = thickness
        self.spar_thickness = thickness # TODO make spar thickness independent of skin-thickness
        self.G_skin = 27 * 10 ** 9  # [Pa] shear modulus
        self.chord = chord
        self.Load = Load
        self.x, self.y = self.loadpoints() * chord
        self.spars = self.createspars(spars)
        # self.cells = self.createcells()
        # self.qT, self.dtheta_dz_T = self.calc_qt()
        self.n = n
        self.indx_lastspar = int(self.spars[-1][-1]) + 1
        self.skin_geom = LineString(np.column_stack((self.x[-self.indx_lastspar:self.indx_lastspar],
                                                     self.y[-self.indx_lastspar:self.indx_lastspar])))
        self.centroid = self.calculate_centriod()
        self.calculate_I()
        # self.Ixx = None
        # self.Iyy = None
        # self.Ixy = None

    def loadpoints(self):
        """ loads airfoil coordinates from file"""
        points = np.array(importAirfoil(self.file))
        plt.plot(points[0], points[1])
        plt.scatter(points[0], points[1])
        plt.axis('equal')
        plt.show()
        return points

    def createspars(self, spars):
        """"
        finding the nodes closesst to the location of the spars
        and linear interpolate to get the ymin and ymax of the spars
        """
        spargeom = np.zeros((len(spars), 4))  # [x, ymin, ymax, indx1]
        for i, x in enumerate(spars * self.chord):
            indx1 = np.searchsorted(self.x[len(self.x) // 2:], x) + len(self.x) // 2
            indx2 = indx1 - 1
            indx3 = -indx1
            indx4 = indx3 - 1
            x1, y1 = self.x[indx1], self.y[indx1]
            x2, y2 = self.x[indx2], self.y[indx2]
            x3, y3 = self.x[indx3], self.y[indx3]
            x4, y4 = self.x[indx4], self.y[indx4]
            ymin = (y2 - y1) / (x2 - x1) * (x - x1) + y1
            ymax = (y4 - y3) / (x4 - x3) * (x - x3) + y3
            spargeom[i] = [x, ymin, ymax, indx1]

            # shifting the points on the airfoil to match spar position
            self.x[indx1] = x
            self.x[indx4] = x
            self.y[indx1] = ymin
            self.y[indx4] = ymax
            ymin = (y2 - y1) / (x2 - x1) * (x - x1) + y1
            ymax = (y4 - y3) / (x4 - x3) * (x - x3) + y3
            spargeom[i] = [x, ymin, ymax, indx1]
            print(indx1)

        return spargeom

    def createcells(self):
        cells = {}
        for i, x in enumerate(self.spars):
            indx = int(x[3])
            cellcount = 1
            if i == 0:
                #first cell
                dsection_y = self.y[-indx - 1:indx + 1]
                dsection_x = self.x[-indx - 1:indx + 1]
                skin = np.column_stack((dsection_x, dsection_y))
                cells[f'cell_{cellcount}'] = cell([skin], [x], [self.thickness], [self.thickness]) # TODO add diffrent thicknesses
                cellcount += 1
                plt.plot(dsection_x, dsection_y, color='pink')
                plt.plot([x[0], x[0]], [x[1], x[2]], color='blue')
                plt.show()

            if i < len(self.spars) - 1:
                #middle cells
                spar2 = self.spars[i + 1]
                indx2 = int(spar2[3])
                bottom = np.column_stack((self.x[indx:indx2 + 1], self.y[indx:indx2 + 1]))
                top = np.column_stack((self.x[-indx2 - 1:-indx], self.y[-indx2 - 1:-indx]))
                cells[f'cell_{cellcount}'] = cell([top, bottom], [x, spar2], [self.thickness, self.thickness], [self.thickness, self.thickness])
                cellcount += 1
                plt.plot(top[:, 0], top[:, 1], color='pink')
                plt.plot([x[0], x[0]], [x[1], x[2]], color='blue')
                plt.plot(bottom[:, 0], bottom[:, 1], color='pink')
                plt.plot([spar2[0], spar2[0]], [spar2[1], spar2[2]], color='red')
                plt.show()
        for i in cells.values():
            i.make_sections(self.sparresolution)

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

        for i, region in enumerate(self.cells.values()):
            if i == 0:
                A[0, 0] = (1 / 2 / region.area / self.thickness / self.G_skin * region.skinperimeter
                           + 1 / 2 / region.area / self.thickness / self.G_skin * region.spar_perimeter[
                               -1])  # TODO add G_spar and tspar as parameters
                A[0, 1] = -1 / 2 / region.area / self.thickness / self.G_skin * region.spar_perimeter[
                    -1]  # TODO add G_spar and tspar as parameters
            else:
                A[i, i] = (
                            1 / 2 / region.area / self.thickness / self.G_skin * region.skinperimeter  # TODO add G_skin_top/bottom and tskin_top/bottom as parameters
                            + 1 / 2 / region.area / self.thickness / self.G_skin * region.spar_perimeter[
                                0]  # TODO add G_spar1 and tspar1 as parameters
                            + 1 / 2 / region.area / self.thickness / self.G_skin * region.spar_perimeter[
                                -1])  # TODO add G_spar2 and tspar2 as parameters
                A[i, i + 1] = -1 / 2 / region.area / self.thickness / self.G_skin * region.spar_perimeter[
                    -1]  # TODO add G_spar2 and tspar2 as parameters
                A[i, i - 1] = -1 / 2 / region.area / self.thickness / self.G_skin * region.spar_perimeter[
                    0]  # TODO add G_spar1 and tspar1 as parameters

            A[i, -1] += -1
            A[-1, i] += 2 * region.area

        qt = np.linalg.solve(A, x)  # last element is dtheta/dz
        return qt[:-1], qt[-1]

    def create_nodes(self):
        pass
    def calculate_I(self):
        """calculate Ixx, Iyy, Ixy, of the section
        """
        self.Ixx = 0
        self.Iyy = 0
        self.Ixy = 0
        xs = self.x[-self.indx_lastspar:self.indx_lastspar] - self.centroid[0]
        ys = self.y[-self.indx_lastspar:self.indx_lastspar] - self.centroid[1]
        for i, _ in enumerate(xs):
            x1 = xs[i]
            x2 = xs[i+1]
            y1 = ys[i]
            y2 = ys[i+1]
            yc = (ys[i] + ys[i+1]) / 2
            xc = (xs[i] + xs[i+1]) / 2
            ds = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
            dA = ds * self.thickness
            yoverx = (y1 - y2) / (x1 - x2)
            beta = np.arctan(yoverx) * 180 / np.pi
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

        pass
    def calculate_centriod(self):
        """
        calculate the centroid of the structural section bounded by the last spar
        :return:
        """
        skincentroid = np.array(self.skin_geom.centroid.coords.xy)
        skinarea = self.skin_geom.length * self.thickness
        spar_centroids =np.zeros((len(self.spars), 3))
        for i, s in enumerate(self.spars):
            spar_lenght = s[2] - s[1]
            spar_y = s[1] + spar_lenght / 2
            spar_area = spar_lenght * self.spar_thickness
            spar_centroids[i] = [s[0], spar_y, spar_area]
        centroid = (skincentroid * skinarea + np.sum(spar_centroids[:, 0:2] * spar_centroids[:, 2])) / (
                skinarea + np.sum(spar_centroids[:, 2]))

        return centroid



class cell():
    def __init__(self, skin_sections: list, spars: list, skin_thickness: float, spar_thickness: list):

        self.skin_sections = skin_sections
        self.spars = spars
        self.skin_thickness = skin_thickness
        self.spar_thickness = spar_thickness
        self.polygon, self.ring = self.make_polygon()
        self.perimeter = self.ring.length
        self.area = self.polygon.area
        self.spar_perimeter = [i[2] - i[1] for i in self.spars]
        self.skinperimeter = self.perimeter - np.sum(self.spar_perimeter)
        self.sections = None
        self.nodes = None

        # print('cell created with area', self.area)
        # print('Total perimeter', self.perimeter)
        # print('Spar perimeter', self.spar_perimeter)
        # print('Skin perimeter', self.skinperimeter)


    def make_polygon(self):
        print(f'amount of sections {len(self.skin_sections)}')
        if len(self.skin_sections) == 1:
            return Polygon(*self.skin_sections), LinearRing(*self.skin_sections)
        else:
            return Polygon(np.stack(self.skin_sections, axis=0).reshape((-1, 2))),LinearRing(np.stack(self.skin_sections, axis=0).reshape((-1, 2)))


    def make_sections(self, spar_resolution):
        self.sections = []
        for i, coord1 in enumerate(self.skin_sections[0]):
            if i != len(self.skin_sections[0]) - 1:
                coord2 = self.skin_sections[0][i+1]
                self.sections.append(section(coord1[0], coord1[1], coord2[0], coord2[1], self.skin_thickness))

        x = self.spars[0][0]
        lenght = self.spars[0][2] - self.spars[0][1]
        n = lenght // spar_resolution + 1
        ys = np.linspace(self.spars[0][1], self.spars[0][2], n, endpoint=True)
        for j,y in enumerate(ys):
            if j != len(ys) - 1:
                self.sections.append(section(x, y, x, ys[j+1], self.spar_thickness[0]))

        if len(self.spars) > 0:
            #sections bottom skin
            for k, coord1 in self.skin_sections[-1]:
                if k != len(self.skin_sections[-1]) - 1:
                    coord2 = self.skin_sections[-1][k + 1]
                    self.sections.append(section(coord1[0], coord1[1], coord2[0], coord2[1], self.skin_thickness))
            #sections right Spar
            x = self.spars[1][0]
            n = lenght // spar_resolution + 1
            ys = np.linspace(self.spars[0][1], self.spars[0][2], n, endpoint=True)
            for l, y in enumerate(ys):
                if l != len(ys) - 1:
                    self.sections.append(section(x, y, x, ys[l + 1], self.spar_thickness[0]))
    def make_nodes(self):
        """
        loop through sections
        calculate ds
        create nodes
        calculate centroid
        add xprime and yprime
        compute Ixx, Iyy, Ixy total cell


        :return:
        """
        for i, section in enumerate(self.sections):
            if i < len(self.sections) - 1:
                area = 0.5 * section.Aeq + 0.5 * self.sections[i + 1].Aeq
                x, y = section.x2, section.y2
            else:
                area = 0.5 * section.Aeq + 0.5 * self.sections[0].Aeq
                x, y = self.section[0].x1, self.section.y1
        pass

class node():
    def __init__(self, x, y, xprime, yprime, A):
        self.x = x # x location in airfoil coordinates
        self.y = y # y location in airfoil coordinates
        self.A = A # Area of the node
        self.Ixx = self.A * yprime**2
        self.Iyy = self.A * xprime**2
        self.Ixy = self.A * xprime * yprime
        self.dqs =  None# jump in shear flow over this node
        self.qs  = None # shear flow due to shear forces on the section ending at this node ccw
    def compute_dqs(self, Ixx, Iyy, Ixy, Vx, Vy):
        pass
class section():
    def __init__(self, x1, y1, x2, y2, thickness):
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.s = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        self.Aeq = self.s * thickness # We assume 0 thickness skin after idealisation


if __name__ == "__main__":
    airfoil = Airfoil("waspairfoil.txt", 0.1, [0.3, 0.5], 1, [0.1, 1, 1])
