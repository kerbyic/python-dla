import numpy as np
import matplotlib.pyplot as plt
from colour import Color
import time

#***DON'T FORGET ACKNOWLEDGEMENTS***#

#*****START PARTICLE CLASS*****#

class particle:
    def __init__(self, coord=(0.0,0.0), radius=1, color=[0,0,1], stickiness = 1, stuck=False):
        self.coord = coord
        self.radius = radius
        self.stuck = stuck
        self.color = color
        self.stickiness = stickiness
    
    # wrappers
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
    
    # add the distance in 2D from add to current coordinates
    def periodic_boundary(self, bnd):
        if self.coord[0] > bnd:
            self.coord[0] -= 2*bnd
        elif self.coord[0] < -bnd:
            self.coord[0] += 2*bnd
        if self.coord[1] > bnd:
            self.coord[1] -= 2*bnd
        elif self.coord[1] < -bnd:
            self.coord[1] += 2*bnd
            
    def move(self, add, bnd):
        self.coord[0] += add[0]
        self.coord[1] += add[1]
        self.periodic_boundary(bnd)
            
    def move_radial(self, deg, step, bnd):
        self.coord[0] += step*np.cos(deg)
        self.coord[1] += step*np.sin(deg)
        self.periodic_boundary(bnd)
            
    # probability to stick the particle based on stickiness
    def stick(self):
        if np.random.rand() < self.stickiness:
            self.stuck = True
        return self.stuck
            
#*****END PARTICLE CLASS*****#

# find if two particle objects are touching
def touching(pnt1, pnt2):
    return (pnt1.get_coord()[0]-pnt2.get_coord()[0])**2 + \
            (pnt1.get_coord()[1]-pnt2.get_coord()[1])**2 < (pnt1.radius + pnt2.radius)**2

# find angle to origin from point
def origin_angle(pnt):
    pnt = -np.asarray(pnt)
    mag = np.sqrt(pnt[0]**2+pnt[1]**2)
    dot = np.dot([1,0],pnt)
    deg = np.arccos(dot/mag)
    if pnt[1] < 0:
        deg += 2*(np.pi-deg)
    return deg

# return any random angle over all 2D space
def any_angle():
    return np.random.rand()*2*np.pi

# return random values according to normal distribution with ability to restrict movement
def norm_variation(scale=1, backward=True):
    var = np.random.normal(scale=scale)
    if not(backward):
        while var < -np.pi/2 or var > np.pi/2:
            var = np.random.normal(scale=scale)
    return var

# return random point along a circle of radius rad
def rplace_circular(rad):
    deg = np.random.rand()*2*np.pi
    return rad*np.array((np.cos(deg), np.sin(deg)))

# return random point along side of square, side half length bound, centered at origin
def rplace_square(bound):
    coord = np.random.rand(2)*2*bound-bound
    coord[np.random.randint(2)] = np.random.choice([-bound,bound])
    return coord

# return random point on a horizontal line going across domain
def rplace_hline(bound, y=0):
    return np.array([np.random.rand()*2*bound-bound, y])

# return random point on a vertical line going up domain
def rplace_vline(bound, x=0):
    return np.array([x, np.random.rand()*2*bound-bound])

def rplace_random(bound):
    return np.random.rand(2)*2*bound-bound

def rplace_semicircle(rad, lower, upper):
    deg = np.random.rand()*(upper-lower)+lower
    return rad*np.array([np.cos(deg), np.sin(deg)])

def rplace_hexagonal(bound):
    deg = np.random.choice(np.linspace(0, 5/3*np.pi, 6))
    return bound*np.array([np.cos(deg), np.sin(deg)])


# find the particle of largest radius
def get_largest_radius(particlelist):
    largest = 0
    for i in particlelist:
        largest = i.radius if i.radius > largest else largest
    return largest

# get the cell index of current position
def get_cell(coord, cell_length, bound):
    x, y = coord[0]+bound, coord[1]+bound
    return [int(x/len_cell), int(y/len_cell)]

# generate seeds in a circle according to the boundary with a specified radius
def seed_circle(bound, radius=1, num=360):
    seeds = []
    for i in range(num):
        x, y = bound*np.cos(i/num*2*np.pi), bound*np.sin(i/num*2*np.pi)
        seeds.append(particle((x,y), color=[1,0,0], stuck=True, radius=radius))
    return seeds

# generate seeds in a vertical line at a given x
def seed_vline(bound, num, x=0, radius=1):
    y_range = np.linspace(-bound,bound,num)
    seeds = []
    for i in y_range:
        seeds.append(particle((x,i), color=[1,0,0], stuck=True, radius=radius))
    return seeds

# generate seeds in a horizontal line at a given y
def seed_hline(bound, num, y=0, radius=1):
    x_range = np.linspace(-bound,bound,num)
    seeds = []
    for i in x_range:
        seeds.append(particle((i,y), color=[1,0,0], stuck=True, radius=radius))
    return seeds

# generate seeds randomly over entire region
def seed_random(bound, num, radius=1):
    seeds = []
    for i in range(num):
        coord = np.random.rand(2)*2*bound-bound
        seeds.append(particle(coord, color=[1,0,0], stuck=True, radius=radius))
    return seeds
    
seed = int(np.random.rand()*1e9)
#seed = 216944555
#seed = 808651310
seed= 259058252
print(seed)
np.random.seed(seed)

# simulation constants
radius = 1
stepsize = 0.1
waves = 30
particles = 600
bound = 40
N = 20*bound/stepsize

# color generation
blue = Color('blue')
gradient = list(blue.range_to(Color('red'), waves))

# generate seeds
cluster = []
cluster.append(particle([0,0], color=[1,0,0], stuck=True, radius=radius))
#cluster = seed_circle(bound, radius, 360)
#cluster = seed_hline(bound, 100, -bound)
#cluster = seed_random(bound, 1)

# generate cell parametrs and create array
# cell length is based on largest in cluster to reduce/eliminate chance
# of missing neighbors that are actually in range
len_cell = max([2*get_largest_radius(cluster)/np.sqrt(2), 2*radius/np.sqrt(2)])*2
ncells = int(2*bound/len_cell)+1
cell_array = np.empty((ncells,ncells), dtype=list)

for i in range(ncells):
    for j in range(ncells):
        cell_array[i,j] = []

# add seeds to correpsonding cell
for i in cluster:
    current_cell = get_cell(i.get_coord(), len_cell, bound)
    cell_array[current_cell[0], current_cell[1]].append(i)

# start simulation
start = time.time()

# loop through waves
for l in range(waves):
    
    # initialize the number of steps for particles before continuing
    steps = N
#    num = 5+l
    free_particle = []
    
    # place particles at coordinates according to desired boundary and append to free_particle
    for k in range(particles):
#        coord = rplace_circular(0.1)
#        coord = rplace_hline(bound, bound)
#        coord = rplace_random(bound)
        coord = rplace_hexagonal(bound)
#        size = radius/(l+1)**(1/3)
        size = radius
        free_particle.append(particle(coord, stuck=False, radius=size, color=gradient[l].rgb, stickiness=1))
    
    # set the all_stuck condition to run loop at least once
    all_stuck = False
    
    # walk the particles, looping through each individually
    while not(all_stuck) and steps > 0:
        # deep copy of the cluster to append to
        temp_cluster = cluster[:]
        # a list of the stuck condition for each particle in the wave
        stuck_list = []
        # loop through particles in wave updating movement and stuck condition
        # do nothing to stuck particles
        for i in free_particle:
            # if the particle is free, move
            if i.stuck == False:
                # movement according to function
                bias = origin_angle(i.get_coord())
#                bias = 1.5*np.pi
#                degree = any_angle()
                degree = norm_variation(3, True) + bias
                i.move_radial(degree, stepsize, bound)
                # get current cell and populate list of nearby seeds to stick to
                current_cell = get_cell(i.get_coord(), len_cell, bound)
                surround_list = []
                for xcell in [0,-1,1]:
                    for ycell in [0,-1,1]:
                        x = current_cell[0]+xcell
                        if x<0:
                            continue
                        elif x>ncells-1:
                            continue

                        y = current_cell[1]+ycell
                        if y<0:
                            continue
                        elif y>ncells-1:
                            continue
                        surround_list.extend(cell_array[x,y])
                # loop through each neighbor and determine distance and stick if touching
                # break once the particle is stuck
                for j in surround_list:
                    if touching(i, j):
                        i.stick()
                        if i.get_stuck():
                            temp_cluster.append(i)
                            cell_array[current_cell[0], current_cell[1]].append(i)
                        break
            
            # update the stuck list for each particle and determine if all are true
            stuck_list.append(i.stuck)
        all_stuck = np.isin(stuck_list, True).all()
        # save the temporary list to the actual cluster list with deep copy
        cluster = temp_cluster[:]
        steps -= 1
    
    # print progress info
    print(len(cluster))
    print("{:.2f}%".format((l+1)/waves*100))

# print cluster and runtime information
print(len(cluster))
print("{:.2f}ms".format((time.time()-start)*1000))

# plot all particles in the cluster
fig = plt.figure()
ax = fig.add_subplot(111)
for i in cluster:
    ax.add_artist(i.get_artist())

ax.set_ylim(-bound, bound)
ax.set_xlim(-bound, bound)
fig.gca().set_aspect('equal', adjustable='box')
