import numpy as np
from matplotlib import pyplot as plt
import single as s 
reload(s)


#plots stable h against D
critical_h1 = []
critical_h2 = []
critical_h3 = []
D = np.arange(0., 1., 0.01)
for i in D:
    pend = s.single(0.1, i) 
    pend1 = s.single((3./4.)*np.pi, i)
    a = pend.stability(pend.rk4, 0.01, 20)     
    b = pend.stability(pend.expl_euler, 0.01, 20)    
    c = pend.stability(pend1.big_angle_euler, 0.01, 20)               
    critical_h1.append(a)
    critical_h2.append(b)
    critical_h3.append(c)
plt.plot(D, critical_h1, label = "RK4")
plt.plot(D, critical_h2, label = "Explicit Euler \n(small angle)")
plt.plot(D, critical_h3, label = "Explicit Euler \n(arbitrary angle)")
plt.xlabel("Scaled damping coefficient $\^D$ [D/(m$\sqrt{gl}$)]")
plt.ylabel("Critical Step Length $\Delta$$\^t$   [t$\sqrt{l/g}$]")
plt.title("Stability Condition and Damping Coefficient")
plt.legend(loc=6, prop={'size':10})
plt.show()       

pend = s.single(0.1, 0.)
pend1 = s.single(0.1, 1.)
pend.plot(pend.expl_euler, 0.04, 100, 1)
pend.plot(pend.leapfrog, 0.04, 100, 2)
pend.plot(pend.impl_euler, 0.04, 100, 3)
pend.plot(pend.rk4, 0.04, 100, 4)

stable_h = pend.stability(pend.expl_euler, 0.0001, 200)
pend.stability(pend.rk4, 2., 200)
pend.plot_energy(pend.expl_euler, 2, 100, 0)
pend.plot_energy(pend.leapfrog, 2., 100, 1)
for i in np.arange(1., 2., 0.1):
    pend.plot_energy(pend.impl_euler, i, 100, 2)
    pend.plot_energy(pend.rk4, i, 100, 3)
pend1.plot_energy(pend1.expl_euler, 1., 100, 0)
pend1.plot_energy(pend1.leapfrog, 2., 100, 1)
pend1.plot_energy(pend1.impl_euler, i, 100, 2)
pend1.plot_energy(pend1.rk4, i, 100, 3)

pend.stability(pend.leapfrog, 0.1, 200)
pend.stability(pend.leapfrog, 1., 200)
for i in np.arange(0.1, 0.5, 0.1):
    pend.plot_energy(pend.big_angle_euler, i, 200)
    pend.plot_energy(pend.expl_euler, i, 200)

