from matplotlib import pyplot as plt
import numpy as np
import double_pendulum as dp
reload(dp)

"""critical = []
critical2 = []
R = np.arange(0.01, 10, 0.1)
for i in R:
    pendulum = dp.double(i, 0.)
    pendulum2 = dp.double(i, 0.2)
    critical.append(pendulum.stability(0.1, 50))
    critical2.append(pendulum2.stability(0.1, 50))

plt.figure()
plt.plot(R, critical, label = "$\^D$ = 0.")
plt.plot(R, critical2, label = "$\^D$ = 0.2")
plt.xlabel("R")
plt.ylabel("Critical Step Length $\Delta$$\^t$   [t$\sqrt{l/g}$]")
plt.legend()
plt.show()"""
    
pend1 = dp.double(0.01, 1.)
pend2 = dp.double(1., 1.)
pend3 = dp.double(100., 1.)
pend4 = dp.double(0.01, 0.)
pend5 = dp.double(1., 0.)
pend6 = dp.double(100., 0.)
pend1.plot(0.01, 100, 1)
pend2.plot(0.01, 100, 2)
pend3.plot(0.01, 20, 3)
#pend1.plot_energy(1., 50)
#pend2.plot_energy(1., 50)
#pend3.plot_energy(1., 50)
#pend4.plot_energy(1., 50)
#pend5.plot_energy(1., 50)
#pend6.plot_energy(0.19, 50)

#pend6.stability(, 50)"""
