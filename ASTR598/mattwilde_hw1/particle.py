class Position:
    
    def __init__(self, xx=0.0, yy=0.0, zz=0.0):
        """Positon measured in meters."""
        
        self.x = xx
        self.y = yy
        self.z = zz

class Velocity:
    
    def __init__(self, vxx=0.0, vyy=0.0, vzz=0.0):
        """velocity in m/s"""
        self.vx = vxx
        self.vy = vyy
        self.vz = vzz
        
class Particle:
    
    def __init__(self, m, Position, Velocity):
        self.mass = m
        self.pos = Position
        self.vel = Velocity
        self.x = Position.x
        self.y = Position.y
        self.z = Position.z
    
    @staticmethod
    def force(p1, p2):
        G = 6.67408e-11
        r2 = ((p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2)
        return G*(p1.mass * p2.mass) / r2
        #return (fx**2 + fy**2 + fz**2)**(1./2.)
