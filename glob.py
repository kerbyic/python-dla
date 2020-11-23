import matplotlib.pyplot as plt

class particle:
    def __init__(self, coord=(0.0,0.0), radius=1, color=[0,0,1], stuck=False):
        self.coord = coord
        self.radius = radius
        self.stuck = stuck
        self.color = color
    def set_stuck(self, stuck):
        self.stuck = stuck
    def get_stuck(self):
        return self.stuck
    def set_coord(self, coord):
        self.coord = coord
    def get_coord(self):
        return self.coord
    def set_radius(self, rad):
        self.radius = rad
    def get_radius(self):
        return self.radius
    def set_color(self, color):
        self.color = color
    def get_color(self):
        return self.color
    def get_artist(self):
        return plt.Circle(self.coord, self.radius, color=self.color)
    def move(self, add, bnd):
        self.coord[0] += add[0]
        self.coord[1] += add[1]
        if self.coord[0] > bnd:
            self.coord[0] -= 2*bnd
        elif self.coord[0] < -bnd:
            self.coord[0] += 2*bnd
        if self.coord[1] > bnd:
            self.coord[1] -= 2*bnd
        elif self.coord[1] < -bnd:
            self.coord[1] += 2*bnd