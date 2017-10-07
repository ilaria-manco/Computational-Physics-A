from matplotlib import pyplot as plt
import numpy as np

class double:
    """Class that represent a double pendulum"""
    
    def __init__(self, R, G):
        """Initialise a pendulum with initial angle and damping constant.
        
        Arguments:
        theta_0 -- initial angle 
        D -- damping constant 
        h -- time increment
        t -- time
        """
        self.G = float(G)
        self.R = R
    
    def y(self, h, time):
        t = np.arange(0, time, h)
        #initial conditions theta_0 = 0.1, phi_0 = 0, velocities = 0
        return np.array([np.full(t.size, 0.1),np.full(t.size, 0.) , np.full(t.size, 0.), np.full(t.size, 0.)])
        
    def matrix(self):
        R = self.R
        G = self.G
        
        mat = np.array([[0., 0., 1., 0.], [0., 0., 0., 1.], [-(R+1.), R, -G, 0.], [R+1., -(R+1.), G*(1.-(1./R)), -G/R]])
        return mat
        
    def rk4(self, h, time):
        t = np.arange(0, time, h)
        M = self.matrix()
        y = self.y(h, time)

        for i in range(1, t.size):
            #np.matrix[:,0] to access column
            k1 = np.dot(M,y[:,i-1])
            y1 = y[:,i-1] + (h/2.)*k1
            k2 = np.dot(M, y1)
            y2 = y[:,i-1] + (h/2.)*k2
            k3 = np.dot(M, y2)
            y3 = y[:,i-1] + (h)*k3
            k4 = np.dot(M, y3)
            y[:,i] = y[:,i-1] + h/6. * (k1 + 2*k2 + 2*k3 + k4)
        
        return y[0,:], y[1,:], y[2,:], y[3,:]
    
    def energy(self, h, time):
        R = self.R
        
        theta = self.rk4(h, time)[0]
        phi = self.rk4(h, time)[1]
        omega = self.rk4(h, time)[2]
        v = self.rk4(h, time)[3]
        
        ke = (omega**2)/2. + (R*omega**2)/2. + (R*v**2)/2. + R*omega*v
        pe = (1./2.)*theta**2 + R/2.*(theta**2 + phi**2)
        energy = ke + pe
        return ke, pe, energy
        
    def plot(self, h, time, i):
        t = np.arange(0, time, h)
        
        theta = self.rk4(h, time)[0]
        phi = self.rk4(h, time)[1]
        
        a = 310 + i
        plt.subplot(a)
        plt.plot(t, theta, label = "$\\theta$")
        plt.plot(t, phi, label = "$\phi$")
        plt.xlabel("$\^t$ [t $\sqrt{g/l}$]")
        title = "R ="+ str(self.R)
        plt.title(title)
        
        e = self.energy(h, time)[2]
        #plt.plot(t, e, label = "$\^E$ [E/mgl]")
        plt.legend(loc=2, ncol = 2, prop={'size':10})
        plt.show()
        
    def plot_energy(self, h, time):
        """Plot log(E/E_0) vs time."""

        t = np.arange(0, time, h)
        energy = self.energy(h, time)[2]
        log_energy = np.log(energy/energy[0])
        plt.xlabel("$\^t$ [t $\sqrt{g/l}$]")
        plt.ylabel("log($\^E$/$\^E_0$)")
        plt.plot(t, log_energy, label = "D="+ str(self.G) + "R="+ str(self.R))
        plt.legend(loc=2, prop={'size':10})
        plt.show()
    
    def stability(self, h_min, time):
        """Return the critical time step length at which the method becomes unstable.
        
        Arguments:
        h_min  -- smallest time step to start testing from
        time   -- time interval to run the test for 
        """
        step = np.arange(h_min, 5., 0.0001)
        if self.G == 0.:                
            for n in step:
                energy = self.energy(n, time)[2]
                energy_in = energy[:(len(energy)-1)/2]
                energy_fin = energy[(len(energy)-1)/2 :]
                #unstable if mean of energy increases or decreases by 0.1%
                if (np.mean(energy_fin) - energy[0])/(energy[0]) > 0.001:
                    print n
                    return n
                    
        if self.G > 0.:
            for n in step:
                energy = self.energy(n, time)[2]

                energy_in = energy[:(len(energy)-1)/2]
                energy_fin = energy[(len(energy)-1)/2 :]
                if np.mean(energy_fin) > np.mean(energy_in):
                    print n
                    return n

