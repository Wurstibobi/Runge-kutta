import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#Annetaan vakiot Au / paiva yksikoissa
g = 2.959*10**-4
m = np.array([1.47,9.2*0.0009543,0.0009543*8.3,7*0.0009543])
r0 = [-0.5688,0]
r1 = [16.4,0]
r2 = [-42.9,0]
r3 = [68,0]
v0 = np.array([0,-0.000025])
v1 = np.array([0,(g*m[0]/r1[0])**(1/2)])
v2 = np.array([0,-1.0001*(g*m[0]/-r2[0])**(1/2)])
v3 = np.array([0,0.9999*(g*m[0]/r3[0])**(1/2)])

s = np.array([r0[0],r0[1],r1[0],r1[1],r2[0],r2[1],r3[0],r3[1],v0[0],v0[1],v1[0],v1[1],v2[0],v2[1],v3[0],v3[1]])

#Luodaan funktio, joka laskee kiihtyvyyden ja nopeuden

def funktio(s,t):
    r0 = s[:2]
    r1 = s[2:4]
    r2 = s[4:6]
    r3 = s[6:8]
    v0 = s[8:10]
    v1 = s[10:12]
    v2 = s[12:14]
    v3 = s[14:16]
    r01 = r1-r0
    r02 = r2-r0
    r03 = r3-r0
    r12 = r2-r1
    r13 = r3-r1
    r23 = r3-r2
    r01s = (r01*r01).sum()**1.5
    r02s = (r02*r02).sum()**1.5
    r03s = (r03*r03).sum()**1.5
    r12s = (r12*r12).sum()**1.5
    r13s = (r13*r13).sum()**1.5
    r23s = (r23*r23).sum()**1.5

    a0 = -g*m[1]*-1*r01/r01s -   g*m[2]*-1*r02/r02s-   g*m[3]*-1*r03/r03s
    a1 = -g*m[0]*r01/r01s -   g*m[2]*-1*r12/r12s -  g*m[3]*-1*r13/r13s
    a2 = -g*m[0]*r02/r02s - g*m[1]*r12/r12s - g*m[3]*-1*r23/r23s
    a3 = -g*m[0]*r03/r03s - g*m[1]*r13/r13s -g*m[2]*r23/r23s

    return np.array([v0[0],v0[1],v1[0],v1[1],v2[0],v2[1],v3[0],v3[1],a0[0],a0[1],a1[0],a1[1],a2[0],a2[1],a3[0],a3[1]])

#Nyt teemme rungekutta funktion
#Kun annamme kiihtyvyyden ja nopeuden se laskee kappaleen paikan aika-askeleen jalkeen


def rk4(s, dt, t):
    k1 = funktio(s, t)*dt
    k2 = funktio(s+k1/2.0, t+dt/2.0)*dt
    k3 = funktio(s+k2/2.0, t+dt/2.0)*dt
    k4 = funktio(s+k3, t+dt)*dt
    return s + (k1+2*k2+2*k3+k4)/6.0


p1x = []
p2x = []
p3x = []
p4x = []
p1y = []
p2y = []
p3y = []
p4y = []
t = 0

dt = 1
i = 0

# Tassa on iteraatiomme, jolle annetaan aikavali paivissa.
# Iteraatio paivittaa kappaleitten paikat
# Sitten koodi lahtee vetamaan uuden kierroksen uusilla paikoilla
while t < 365*100:
    s = rk4(s,dt,t)
    i = i+1
    t = t +dt
    if (i%100 == 0):
        p1x.append(s[0])
        p1y.append(s[1])
        p2x.append(s[2]-s[0])
        p2y.append(s[3]-s[1])
        p3x.append(s[4]-s[0])
        p3y.append(s[5]-s[1])
        p4x.append(s[6]-s[0])
        p4y.append(s[7]-s[1])
        print(s[1])


#Lopusi plottaamme kappaleitten paikat

plt.axes().set_aspect('equal')
plt.ylim(-70,70)
plt.xlim(-70,70)
plt.title('1')
plt.xlabel('X(Au)')
plt.ylabel('Y(Au)')
plt.plot(r1[0],r1[1],'r.',markersize=10)
plt.plot(r2[0],r2[1],'c.',markersize=10)
plt.plot(r3[0],r3[1],'k.',markersize=10)
plt.plot(0,0, 'b.',label = 'Aurinko',markersize = 10)
plt.plot(p2x,p2y, 'r.',label = 'Planeetta E',markersize = 3)
plt.plot(p3x,p3y, 'c.',label = 'Planeetta C',markersize = 3)
plt.plot(p4x,p4y, 'k.',label = 'Planeetta B',markersize = 3)
plt.legend()
plt.show()
