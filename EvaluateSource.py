from MakeLightcurves import *

print(get_star_name(tic=69747919))

tess_lc = make_lc(tic=69747919,mission='tess')
k2_lc = make_lc(tic=69747919,mission='k2')

plot_lc(tess_lc,tic=69747919)
plot_lc(k2_lc,tic=69747919)
