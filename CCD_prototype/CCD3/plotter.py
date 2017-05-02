#this is a simple plotter where you manually insert data sets (in case you are interested in spesific cases)

import itertools
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc

xData = range(2,42)
yData = [14.474197721372551,
		 14.218743793305565,
		 14.207942941796345,
		 13.936112367402592,
		 13.611507926048816,
		 13.611507926048816,
		 13.57451376342269,
		 13.534148280969257,
		 13.505725917938069,
		 13.494319125497038,
		 13.493846612339542,
		 13.492441495335218,
		 13.49155729440426,
		 13.49155729440426,
		 13.49146272954267,
		 13.491331913824345,
		 13.491295036250108,
		 13.491289243495189,
		 13.491285320314288,
		 13.491277400922622,
		 13.49127102905164,
		 13.49127102905164,
		 13.491267741663021,
		 13.491265278557359,
		 13.49126066855264,
		 13.491259278823192,
		 13.491259278823192,
		 13.491257416080982,
		 13.491256539828935,
		 13.491256539828935,
		 13.491256450040053,
		 13.491256192889837,
		 13.491256035877987,
		 13.491255921457979,
		 13.491255880376345,
		 13.491255847030002,
		 13.491255798227725,
		 13.491255798227725,
		 13.491255789698593,
		 13.491255774417763]
print len(xData)


#different marker for each graph
marker = itertools.cycle(('s', 'v', 'o', 'h'))
linestyle = itertools.cycle(('-', '--', '-.', ':'))

#plot commands
fig, ax = plt.subplots()

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax.grid(True)

#iterate over all data sets
#counter = 0
#for ECC in Vals:
#	plt.plot(Rho, ECC, label=r"$N_h = {}$".format(Nh[counter]), marker=marker.next(), linestyle = linestyle.next())
#	plt.hold('on')
#	counter += 1

plt.plot(xData, yData, label=r"$N_h = 14$", marker=marker.next(), linestyle = linestyle.next())

plt.legend(loc=4)

plt.xlabel(r"$\rho$ [fm$^{-3}$]", fontsize=16)
plt.ylabel(r"$E_{CCD}/A$ [MeV]", fontsize=16)

plt.xlim(0, 42)

plt.savefig('states_for_MP.pdf', format='pdf')

plt.show()
plt.tight_layout()
