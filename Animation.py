from matplotlib import pyplot as plt
import numpy as np
import matplotlib.animation as animation

R = 1.
G = 0.

h = 0.1
t = np.arange(0, 20, h)

y = np.array([np.full(t.size, 1.5),np.full(t.size, 1.) , np.full(t.size, 0.), np.full(t.size, 0.)])
M = np.array([[0., 0., 1., 0.], [0., 0., 0., 1.], [-(R+1.), R, -G, 0.], [R+1., -(R+1.), G*(1.-(1./R)), -G/R]])

theta = y[0,:]
phi = y[1,:]
omega = y[2,:]
v = y[3,:]

for i in range(1, t.size):
    #np.matrix[:,0] to access column
    k1 = np.dot(M,y[:,i-1])
    y1 = y[:,i-1] + (h/2.)*k1
    k2 = np.dot(M, y1)
    y2 = y[:,i-1] + (h/2.)*k2
    k3 = np.dot(M, y2)
    y3 = y[:,i-1] + (h/2.)*k3
    k4 = np.dot(M, y3)
    y[:,i] = y[:,i-1] + h/6. * (k1 + 2*k2 + 2*k3 + k4)


"""s = np.sin(theta) + np.sin(phi)
p = -np.cos(theta) - np.cos(phi)
plt.subplot(212)
plt.plot(s, p, '-')"""


#code for animation from matplotlib.org

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2.5, 2.5), ylim=(-2.5, 2.5))
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text


def animate(i):
    thisx = [0, np.sin(theta[i]), np.sin(phi[i]) + np.sin(theta[i])]
    thisy = [0, -np.cos(theta[i]), -np.cos(phi[i]) - np.cos(theta[i])]

    line.set_data(thisx, thisy)
    time_text.set_text(time_template % (i*h))
    return line, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(theta)),
                              interval=5, blit=True, init_func=init)

plt.show()