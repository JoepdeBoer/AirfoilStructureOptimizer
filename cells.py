import numpy as np
from shapely import LinearRing, Polygon
class Cell():
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


    # def make_sections(self, spar_resolution):
    #     self.sections = []
    #     for i, coord1 in enumerate(self.skin_sections[0]):
    #         if i != len(self.skin_sections[0]) - 1:
    #             coord2 = self.skin_sections[0][i+1]
    #             self.sections.append(section(coord1[0], coord1[1], coord2[0], coord2[1], self.skin_thickness))
    #
    #     x = self.spars[0][0]
    #     lenght = self.spars[0][2] - self.spars[0][1]
    #     n = lenght // spar_resolution + 1
    #     ys = np.linspace(self.spars[0][1], self.spars[0][2], n, endpoint=True)
    #     for j,y in enumerate(ys):
    #         if j != len(ys) - 1:
    #             self.sections.append(section(x, y, x, ys[j+1], self.spar_thickness[0]))
    #
    #     if len(self.spars) > 0:
    #         #sections bottom skin
    #         for k, coord1 in self.skin_sections[-1]:
    #             if k != len(self.skin_sections[-1]) - 1:
    #                 coord2 = self.skin_sections[-1][k + 1]
    #                 self.sections.append(section(coord1[0], coord1[1], coord2[0], coord2[1], self.skin_thickness))
    #         #sections right Spar
    #         x = self.spars[1][0]
    #         n = lenght // spar_resolution + 1
    #         ys = np.linspace(self.spars[0][1], self.spars[0][2], n, endpoint=True)
    #         for l, y in enumerate(ys):
    #             if l != len(ys) - 1:
    #                 self.sections.append(section(x, y, x, ys[l + 1], self.spar_thickness[0]))
    # def make_nodes(self, Mx, Ixx, Iyy, Ixy):
    #     """
    #     loop through sections
    #     calculate ds
    #     create nodes
    #     calculate centroid
    #     add xprime and yprime
    #     compute Ixx, Iyy, Ixy total cell
    #
    #
    #     :return:
    #     """
    #     for i, section in enumerate(self.sections):
    #         if i < len(self.sections) - 1:
    #             # sig = MxIyy * y - Mxxy * x/ IyyIxx - Ixy^2
    #             # 1/2 = MxIyy * y1 - MxIxy * x1  / MxIyy * y2 - MxIxy * x2
    #             sig1_sig2 =
    #             area += section.Aeq/6 *
    #             x, y = section.x2, section.y2
    #         else:
    #             area = 0.5 * section.Aeq + 0.5 * self.sections[0].Aeq
    #             x, y = self.section[0].x1, self.section.y1
    #     pass