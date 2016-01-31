
import particle


p1 = particle.Particle(2, particle.Position(0,0,0), particle.Velocity(2,2,2))
p2 = particle.Particle(1, particle.Position(1,1,1), particle.Velocity(2,2,2))
f = particle.Particle.force(p1, p2)

print "The force between particle1 and particle2 is {:.2} Newtons".format(f)
