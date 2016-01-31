"""
author: Matt Wilde
email: mwilde@uw.edu

ASTR507 HW#1 particle in a box

"""
# import timing
from numpy import sin, cos, pi, sqrt
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.animation as animation
# from numba import double
from numba.decorators import jit
from scipy.stats import chisquare
# from math import exp

pl.ioff()
# runtime conditionals
run = True
makeMovie = True
physical_radius = True
make_hist = True
twomasses = False
onemass = True
one_b = True
three_c = False
two_c = False
gravity = False
ndensity = False
crosssection = True

# output prefix
pre = None


# initial values runtime values
numParticles = 100
numTsteps = 300
width = 11.0
height = 11.0
radius = 0.02*width  # 2% the size of the box
v0 = 0.1*width
mass = 1.0
if twomasses:
    mass1 = 1.0
    mass2 = 1.0*10.0
else:
    mass1 = mass2 = mass
t_frac = 0.1
time_step = t_frac*(radius/(numParticles*v0))
#time_step = (height*width)/(numParticles*pi*radius**2*v0)
Fg = -9.8*100 # cm/s**2



# define classes for our objects
class Particle():
    def __init__(self, position, radius, mass, velocity):

        self.x = position[0]
        self.y = position[1]
        self.position = np.array([self.x, self.y])
        self.r = radius
        self.m = mass
        self.vx = velocity[0]
        self.vy = velocity[1]
        self.vmag = sqrt(self.vx**2 + self.vy**2)
        self.v = np.array([self.vx, self.vy])
        self.E = 0.5*self.m*(self.vx**2 + self.vy**2)


# define a grid for the particles to move around on
# allow it to change size
class Grid():
    # initialize the grid
    def __init__(self, width, height):
        self.height = height
        self.width = width
        self.particles = []

    # add particles to grid, chosing the number,radius, and mass
    def add_particles(self, num, radius, mass):
        pos = []
        for i in range(num):
            """ pos = [(width-radius)*np.random.random(),(height-radius)*np.
                                                            random.random()]
                pos = [abs(self.width-radius-v0*time_step)*random.random(),abs
                        (self.height-radius-v0*time_step)*random.random()]
            """
            # make a grid of particles so none are overlapping
            for j in range(int(num/10)):
                for k in range(int(num/10)):
                    pos.append([float(j+1.0), float(k+1.0)])
            theta = 2*pi*np.random.random()
            vel = [v0*cos(theta), v0*sin(theta)]
            if twomasses:
                if i % 2 == 0:
                    mass = mass*10.0
                else:
                    mass = 1.0
            else:
                mass = 1.0
            radius = radius
            self.particles.append(Particle(pos[i], radius, mass, vel))


@jit
def distance(p_i, p_j):
    '''A funciton to define the distance between one particle and the
        rest of the particles. Returns distance as a float.
    '''
    return np.sqrt((p_i[0]-p_j[0])**2 + (p_i[1]-p_j[1])**2)


#===================================================
# maxwell boltzman distribution
#===================================================
def MB_v(v, m):
    """Generic maxwell boltzmann distrobution."""
    c = m/E_0
    #c = 6200 #scaling factor
    fp = v*np.exp((-0.5*m*v**2)/(E_0))
    return c*fp


#===================================================
# chi2
#===================================================
def relaxed(speed, m, t):
    """ A function to get difference between measuered dn/dv and
        the value from maxwell-boltzmann. Possibly depricated...
    """
    # threshold = 100
    if t > 0:
        h, bins = np.histogram(speed, bins=5, normed=True)
        center = (bins[:-1] + bins[1:])/2
        #chi2 = ((MB_v(center, m) - h)**2)/(MB_v(center, m))
        chi = chisquare(h,MB_v(center, m))

        return chi


def chi2_scipy_plot():
    "Use built in chi2 funciton from scipy. Return None"""
    chi_list = []        
    for t in range(1,numTsteps-1):
        a = relaxed(v_all[t],mass,t)
        chi_list.append(a[0])
    
    # make plot to show convergence
    t = np.arange(1,numTsteps-1)
    pl.figure(3)
    pl.ylim(0,3.0)
    pl.plot(t[50:numTsteps],chi_list[50:numTsteps])
    s = "N={}, r={}".format(numParticles,radius)
    pl.title(r"$\chi^{2}$") 
    pl.text(500,2,s)
    pl.xlabel("time_step")
    pl.savefig("chi2_N{}_r{}.png".format(numParticles,radius))
    return None
        
    


def Pressure(m, vxf, vxi):
    """Calculate 2D Pressure to check if P = nkT."""
    p = (m/height)*((vxf - vxi)/time_step)
    return p


# make a grid
g = Grid(width, height)

# add particles to grid
g.add_particles(numParticles, radius, mass)


# make empty np.arrays to store all physical data
# np is easier to manipulate than lists
x_0 = np.empty((numParticles), dtype=np.float)
y_0 = np.empty((numParticles), dtype=np.float)
x_t = np.empty((numTsteps, numParticles), dtype=np.float)
y_t = np.empty((numTsteps, numParticles), dtype=np.float)
v_0 = np.empty((numParticles), dtype=np.float)
v_all = np.empty((numTsteps, numParticles), dtype=np.float)
E_all = np.empty((numTsteps, numParticles), dtype=np.float)
m_j = np.empty((numParticles), dtype=np.float)
E_j = np.empty((numParticles), dtype=np.float)
E_t = np.empty((numTsteps, numParticles), dtype=np.float)
flag = np.empty(numParticles, dtype=float)
list_of_chi2 = []
p = np.empty((numTsteps,numParticles), dtype=float)
# fill arrays with initial values at t=0
for j in range(numParticles):
    x_0[j] = g.particles[j].x
    y_0[j] = g.particles[j].y
    v_0[j] = g.particles[j].vmag
    m_j[j] = g.particles[j].m
    E_j[j] = 0.5*(m_j[j]*v_0[j]**2)

# total E
E_0 = E_j.sum()/numParticles 
#if run:
@jit
def evolve():
    """ This is the main loop for running the particles. Returns None """
    dp = 0.    
    # make loop
    for t in range(numTsteps):
        for j in range(numParticles):
            p1 = g.particles[j]
            x1 = p1.position
            v1 = p1.v
            r1 = p1.r
            m1 = p1.m

            # reflect balls off wall
            # track last thing ball hit with flag
            if x1[0]+r1 >= width and flag[j] != numParticles+1:
                flag[j] = numParticles+1
                v1[0] = -v1[0]
                dp += 2*m1*abs(v1[0])
            if x1[0]-r1 <= 0.001 and flag[j] != numParticles+2:
                v1[0] = -v1[0]
                flag[j] = numParticles+2
                dp += 2*m1*abs(v1[0])
            if x1[1]+r1 >= height and flag[j] != numParticles+3:
                if gravity:                
                    v1[1] = v1[1]
                else:
                    v1[1] = - v1[1]
                    flag[j] = numParticles+3
                    dp += 2*m1*abs(v1[1])
            if x1[1]-r1 <= 0.001 and flag[j] != numParticles+4:
                v1[1] = -v1[1]
                flag[j] = numParticles+4
                dp += 2*m1*abs(v1[1])

            x_t[t,j] = g.particles[j].position[0]
            y_t[t,j] = g.particles[j].position[1]

            # E _0[t,j] = g.particles[j].E
            for k in range(j+1, numParticles):
                p2 = g.particles[k]
                v2 = p2.v
                x2 = p2.position
                m2 = p2.m
                r2 = p2.r
                
                # check for collision
                if distance(x1,x2) < r1+r2 and flag[j] != k:
                    # set flag to make sure partilce doesnt get stuck
                    flag[j] = k
                    
                    # mass coefficient
                    M1 = (2*m1)/(m1+m2)
                    M2 = (2*m2)/(m1+m2)

                    # collide particles and get new velocity
                    p1.v = v1-(M1*(np.inner((v1-v2),(x1-x2))/(np.inner((x1-x2),(x1-x2)))))*(x1-x2)
                    p2.v = v2-(M2*(np.inner((v2-v1),(x2-x1))/(np.inner((x2-x1),(x2-x1)))))*(x2-x1)

                # add downward force for gravity
                if gravity:
                    p2.v[1] += 0.5*(Fg)*time_step**2
                    p1.v[1] += 0.5*(Fg)*time_step**2

                # update particle positions
                x2 += p2.v*time_step
                x1 += p1.v*time_step
                
            p[t,j] = dp    
            v_all[t,j] = sqrt(p1.v[0]**2 + p1.v[1]**2)
            E_all[t,j] = 0.5*m1*(p1.v[0]**2 + p1.v[1]**2)
        # check to numerically end loop when equilibrium reached
        list_of_chi2.append(relaxed(v_all[t], m1, t))
        # check pressure        
    E = (numParticles/(width*height)*E_all[numTsteps-1].mean()) 
    P = (numParticles * dp/(numTsteps*width*height))
    print(dp)
    print("p = N*delta_p/delta_t/(w*h) = {0} = nkT(E_avg) = {1}".format(P,E))
    # print("but p = E if p = numParticles*2p = {}".format(dp/(numParticles)))
        # return None
    
#===================================================
# run particle evolution
#===================================================
if run:
    evolve()
    chi2_scipy_plot()
    

if one_b:
    # make four plots with position and velocity distribution
    vzero = (numTsteps*0.01)*np.hstack(v_0)
    # vend = np.hstack(v_all[-int((numTsteps*0.01)):])
    vend = v_all[numTsteps-1]
    pl.clf()
    
    f, ((ax1, ax2), (ax3, ax4)) = pl.subplots(2, 2)
    s = "{} particles after {} timesteps with radius {}".format(numParticles,numTsteps,radius)
    pl.title(s)
    ax1.set_xlim(0., width)
    ax1.set_ylim(0., height)
    ax2.set_xlim(0., width)
    ax2.set_ylim(0., height)
    ax1.scatter(x_0,y_0)
    ax2.scatter(x_t[numTsteps-1],y_t[numTsteps-1])
    ax3.hist(v_0, histtype='step')#, normed=True, color = 'red',label='t=0')
    ax4.hist(vend, histtype='step', normed=True, color = 'black',label='t=end')
        
        
    pl.show()
    if twomasses:
        f.savefig('4plots_{0}_twomasses'.format(numTsteps))
    else:
        f.savefig('4plots_{0}_onemass'.format(numTsteps,))

def mkmovie(x_t,y_t):
    """ A function to make a movie of the particles in a box. Takes 
        the output of evolve() which are the x and y data from the
        simulation. Saves movie to file and returns None.
    """
    fig = pl.figure(1)
    ax = pl.axes(xlim=(0., width), ylim=(0., height))

    if physical_radius:
        @jit
        def init():
            ax.set_xlim(0., width)
            ax.set_ylim(0., height)
            for j in range(numParticles):
                circle = pl.Circle((x_t[0,j],y_t[0,j]), radius=radius)
                ax.add_patch(circle)
        @jit
        def animate(i):
            #print "Animate step " + str(i)
            ax.cla()
            ax.set_xlim(0., width)
            ax.set_ylim(0., height)
            for j in range(numParticles):
                circle = pl.Circle((x_t[i,j], y_t[i,j]), radius=radius)
                ax.add_patch(circle)

    else:
        def init():
            scat = pl.scatter([0.],[0.])
            return scat,
        def animate(i):
            #print "Animate step " + str(i)
            ax.cla()
            ax.set_xlim(0., width)
            ax.set_ylim(0., height)
            scat = pl.scatter(x_t[i],y_t[i])
            return scat,
    
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=numTsteps, interval=100, blit=True)
    if onemass:    
        anim.save('PiB_onemass_{0}_N{1}.mp4'.format(numTsteps,numParticles))
    if twomasses:
        anim.save('PiB_twomasses_{0}_N{1}.mp4'.format(numTsteps,numParticles))
    if gravity:
        anim.save('PiB_with_gravity.mp4')
    return None


if makeMovie:
    mkmovie(x_t,y_t)

def two_b(N,nbins):
    """ 2b) make a historgram averaged over last N timesteps to show convergence
        with MB. Input is N. Returns None. Works best with ~(30,5)."""
    vend_ave = v_all[-N:]
    h,bins = np.histogram(vend_ave, bins=nbins, normed=True)
    center = (bins[:-1] + bins[1:])/2
    lab1 = 'avg over last {} tsteps'.format(N)
    pl.plot(center,h, 'black', label=lab1)

    # plot MB on top
    pl.plot(v,MB_v(v,mass1),'green',label='MB(v,mass1)')
      
    
    # calculate Chi^2
    #chi2 = (((MB_v(center, mass1) - h)**2)/MB_v(center,mass1)).mean()
    a = chisquare(h,MB_v(center,mass))
    
    # makeplot    
    pl.legend()
    pl.ylabel("dN/dv")
    pl.xlabel("v")
    pl.title("2b: data v MB afer {} steps".format(numTsteps)) 
    pl.text(1.5,0.6,"chi2_m1={:0.2}".format(a[0]))
    fname = '2b_one_mass_{}_steps'.format(numTsteps) 
    pl.savefig(fname)
    return None
    
def two_c(N,nbins,mass1,mass2):
    """ 2c) make a historgram averaged over last N timesteps to show convergence
        with MB for both masses. Input is N. Returns None. Works best with ~(30,5)."""
    if twomasses:        
        vend_ave = v_all[-N:]
        h,bins = np.histogram(vend_ave, bins=nbins, normed=True)
        center = (bins[:-1] + bins[1:])/2
        lab1 = 'ave over last {} tsteps'.format(N)
        pl.plot(center,h, 'black', label=lab1)
    
        # plot MB on top
        pl.plot(v,MB_v(v,mass1),'green',label='MB(v,mass1)')
        pl.plot(v,MB_v(v,mass2),'red',label='MB(v,mass2)')
          
        
        # calculate own Chi^2
        #chi21 = (((MB_v(center, mass1) - h)**2)/MB_v(center,mass1)).mean()
        #chi22 = (((MB_v(center, mass2) - h)**2)/MB_v(center,mass2)).mean()
        
        #built in chisquare        
        a = chisquare(h,MB_v(center,mass1))
        b = chisquare(h,MB_v(center,mass2))
        pl.text(1.5,0.6,"chi2_m1={:0.2}\nchi2_m2={:0.2}".format(a[0],b[0]))
        
        # makeplot    
        pl.legend()

        pl.title("2c: two masses data v MB afer {0} steps".format(numTsteps)) 
        fname = '2c_two_masses_{}_steps'.format(numTsteps) 
        pl.savefig(fname)
        return None
    else:
        print("must run with 2 masses")
    
def three_b(N,nbins,mass):
    """ vary the number density and cross section"""
    pass
        
    

if make_hist:
    f = pl.figure(2)
    vstart = np.hstack(v_all[0])
    vend = v_all[numTsteps-1]
    vmid = v_all[int((numTsteps-1)/2)]
    a = f.gca()
    a.set_ylim(0,1)
    s_end = 't={}'.format(numTsteps)
    s_mid = 't={}'.format(int(numTsteps/2))
    pl.hist(vend, histtype='step', normed=True, color = 'black',label=s_end)
    #pl.hist(vmid, histtype='step', normed=True, color = 'green',label=s_mid)
    pl.hist(v_0, histtype='step', normed=False, color = 'red',label='t=1')
    #pl.hist(vzero, histtype='step', normed=False, color = 'm',label='t=initial')
    v = np.arange(0,v_all[numTsteps-1].max(),0.001)
    pl.plot(v,MB_v(v,mass1),'orange',label='MB(v,mass1)')
    if twomasses:
        pl.plot(v,MB_v(v,mass2),'magenta',label='MB(v,mass2)')    
    
    if twomasses:
        pl.title('Two Masses after {} time steps'.format(pre, numTsteps))
        f.savefig('histogram_twomasses_{}_{}.png'.format(numTsteps,t_frac))
    if onemass:
        pl.title('One mass after {} time steps'.format(numTsteps))
        f.savefig('histogram_onemass_{}_{}.png'.format(numTsteps,t_frac))
    if ndensity:
        pl.title('{} particles / volume after {} time steps'.format(numParticles,numTsteps))
        f.savefig('histogram_ndensity_{}_{}.png'.format(numTsteps,t_frac))
    if crosssection:
        pl.title('cross section = {}r0 with {} particles'.format(radius,numParticles))
        f.savefig('histogram_crossection_{}_{}.png'.format(numTsteps,t_frac))
    pl.show()
# calculate initial and final energy
print("E initial:",E_0,"E_final",E_all[numTsteps-1].mean())
print(time_step, numTsteps)
# calculate pressure
P = (numParticles/ (4*width)) * p[numTsteps-1,0]/numTsteps
E = E_all[numTsteps-1].mean()* (numParticles/width**2)
print("P = {}, nkT = N/A * E_ave = {}".format(P,E))



