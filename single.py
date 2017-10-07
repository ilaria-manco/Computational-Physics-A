from matplotlib import pyplot as plt
import numpy as np

class single:
    """Class that represent a single pendulum"""
    
    def __init__(self, theta_0, D):
        """Initialise a pendulum with initial angle and damping coefficient.
        
        Arguments:
        theta_0 -- initial angle 
        D -- damping coefficient 
        h -- time increment
        t -- time
        """
        self.theta_0 = theta_0
        self.D = D
        
    def u(self, h, time):
        """Return the initial angular displacement theta_0 in an array of size t."""
        t = np.arange(0, time, h)
        return np.full(t.size, self.theta_0)
    
    def v(self, h, time):
        """Return the initial angular velocity in an array of size t."""
        t = np.arange(0, time, h)
        return np.zeros(t.size)
            
    def energy(self, method, h, t):
        """Return an array of the energy calculated using one of the FMDs."""
        pe = 1./2.*(method(h, t)[0])**2
        ke = (1./2.)*(method(h, t)[1])**2
        e = pe + ke
        return e
    
    def expl_euler(self, h, time):
        """Return two arrays with the angular displacement and its rate of change
        calculated using the explicit Euler method.""" 
        v = self.v(h, time)
        u = self.u(h, time) 
        D = self.D
        t = np.arange(0, time, h)
        
        for i in range(1, t.size):
            v[i] = v[i-1] - (u[i-1] + D * v[i-1])*h
            u[i] = u[i-1] + v[i-1]*h
        
        return u, v, "Explicit Euler"
    
    def big_angle_euler(self, h, time):
        v = self.v(h, time)
        u = self.u(h, time) 
        D = self.D
        t = np.arange(0, time, h)
        
        for i in range(1, t.size):
            v[i] = v[i-1] - (np.sin(u[i-1]) + D * v[i-1])*h
            u[i] = u[i-1] + v[i-1]*h
        
        return u, v, "Explicit Euler (arbitrary amplitude)"
    
    def leapfrog(self, h, time):
        """Return two arrays with the angular displacement and its rate of change
        calculated using the leapfrog method.""" 
        v = self.v(h, time)
        u = self.u(h, time) 
        D = self.D
        t = np.arange(0, time, h)
        
        for i in range(1, t.size):
    	    v[1] = v[0] - (u[0] + D*v[0])*h
            u[1] = u[0] + v[0]*h
            v[i] = v[i-2] - 2*(u[i-1] + D*v[i-1])*h
            u[i] = u[i-2] + 2*v[i-1]*h
        
        return u, v, "Leapfrog"
    
    def reverse_leapfrog(self, h, time):
        #not fully tested
        """Return two arrays with the angular displacement and its rate of change
        calculated using the leapfrog method with a reversed time direction.""" 
        v = self.v(h, time)
        u = self.u(h, time) 
        D = self.D
        a = np.arange(0, time, h)
        t = a[::-1]
        
        for i in range(1, t.size):
    	    v[1] = v[0] - (u[0] + D*v[0])*h
            u[1] = u[0] + v[0]*h
            v[i] = v[i-2] - (u[i-1] + D*v[i-1])*h
            u[i] = u[i-2] + v[i-1]*h
        
        return u, v, "Reversed Leapfrog"
        
    def impl_euler(self, h, time):
        """Return two arrays with the angular displacement and its rate of change
        calculated using the implicit Euler method.""" 
        v = self.v(h, time)
        u = self.u(h, time) 
        D = self.D
        t = np.arange(0, time, h)
        
        for i in range(1, t.size):
            v[i] = (v[i-1] - u[i-1]*h)/(1 + h*D + h**2)
            u[i] = u[i-1] + v[i]*h
        
        return u, v, "Implicit Euler"
    
    def rk4(self, h, time):
        """Return two arrays with the angular displacement and its rate of change
        calculated using the RK4 method.""" 
        v = self.v(h, time)
        u = self.u(h, time) 
        D = self.D
        t = np.arange(0, time, h)

        for i in range(1, t.size):
            fv1 = - (D*v[i-1] + u[i-1])
            fu1 = v[i-1]
            v1 = v[i-1] + fv1 * h/2.
            u1 = u[i-1] + fu1 * h/2.
            fv2 = -(D*v1 + u1)
            fu2 = v1
            v2 = v[i-1] + fv2 * h/2.
            u2 = u[i-1] + fu2 * h/2.
            fv3 = -(D*v2 + u2)
            fu3 = v2
            v3 = v[i-1] + fv3 * h
            u3 = u[i-1] + fu3 * h
            fv4 = -(D*v3 + u3)
            fu4 = v3
            v[i] = v[i-1] + (fv1 + 2*fv2 + 2*fv3 + fv4)/6. * h
            u[i] = u[i-1] + (fu1 + 2*fu2 + 2*fu3 + fu4)/6. * h  
        
        return u, v, "RK4"
    
    def plot(self, method, h, time, i):
        """Plot of theta and energy calculated using one of the FDMs."""
        t = np.arange(0, time, h)
        
        a = 410 + i
        plt.subplot(a)
        plt.xlabel("$\^t$ [t $\sqrt{g/l}$]")
        plt.plot(t, method(h, time)[0], label="$\\theta$   [rad]")

        energy = self.energy(method, h, time)
        plt.plot(t, energy, label = "Scaled Energy $\^E$ [E/mgl]", color = 'r')
        plt.legend(loc = 1, ncol=2, borderaxespad=0., prop={'size':10})
        plt.title(method(h, time)[2])    
        plt.show() 
        
    def stability(self, method, h_min, time):
        """Return the critical time step length at which the method becomes unstable.
        
        Arguments:
        method -- FDM to test
        h_min  -- smallest time step to start testing from
        time   -- time interval to run the test for 
        """
        step = np.arange(h_min, 10., 0.0001)
        if self.D == 0.:                
            for n in step:
                energy = self.energy(method, n, time)
                energy_in = energy[:(len(energy)-1)/2]
                energy_fin = energy[(len(energy)-1)/2 :]
                if method == self.leapfrog:
                    if np.amax(energy_fin) - np.amin(energy) > 2.*(np.amax(energy_in) - np.amin(energy_in)):
                        print "a", n
                        return n
                else:
                    #unstable if mean of energy increases or decreases by 0.1%
                    if (np.mean(energy_fin) - energy[0])/(energy[0]) > 0.001:
                        print n
                        return n
                    
        if self.D > 0.:
            for n in step:
                energy = self.energy(method, n, time)
                energy_in = energy[:(len(energy)-1)/2]
                energy_fin = energy[(len(energy)-1)/2 :]
                if np.mean(energy_fin) > np.mean(energy_in):
                    print method(n, time)[2], n
                    return n
                    
    def plot_energy(self, method, h, time, color):
        """Plot log(E/E_0) vs time."""
        
        t = np.arange(0, time, h)
        energy = self.energy(method, h, time)
        log_energy = np.log(energy/energy[0])
        plt.xlabel("$\^t$ [t $\sqrt{g/l}$]")
        plt.ylabel("log($\^E$/$\^E_0$)")
        col = ["r", "b", "g", "c"]
        if self.D > 0.:
            plt.plot(t, log_energy, '--', color = col[color])
        else:
            plt.plot(t, log_energy, color = col[color], label=(method(h, time)[2]))
        plt.legend(loc=2)
        plt.show()

def round_error_leapfrog(theta_0, D, h, t):
    #calculates round-off error of leapfrog method
    #not fully tested
    pend = single(theta_0, D)
    #pend.plot(pend.leapfrog, h, t*4*np.pi)
    a = pend.leapfrog(h, 4*np.pi)[0][-1]
    pend = single(a, D)
    #pend.plot(pend.reverse_leapfrog, h, t*4*np.pi)
    b = pend.reverse_leapfrog(h, 4*np.pi)[0][-1]
    #return abs((theta_0 - b)/theta_0)
    return abs(theta_0 - b)
    
