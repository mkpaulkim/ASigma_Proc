import numpy as np
from tkinter import filedialog
import pycubelib.files_functions as ff
import pycubelib.general_functions as gf
import pycubelib.plotting_functions as pf
import pycubelib.data_functions as df

pi = np.pi
pi2 = pi * 2
pi2limit = (-pi, pi)
sxy = (.35, .35)

txt_path = filedialog.askopenfilename(title='TXT file path', filetypes=[('txt files', '*.txt')])
notes, txt_path = ff.read_txt(txt_path)
nx = gf.find_param(notes, 'nx', int)
ny = gf.find_param(notes, 'ny', int)
dx = gf.find_param(notes, 'dx_um', float)
nh = gf.find_param(notes, 'nh', int)
lam_ns = gf.find_param(notes, 'lam_ns', float)
print(f'> txt_path = {txt_path}')
print(f'> notes: \n{notes}')
blank = np.zeros((ny, nx))

hhhp = []
for n in range(nh+1):
    if n == 0:
        hh_path = txt_path.replace('.txt', f'_aa.png')
        hh, _ = ff.read_png(hh_path, alimit=(0., 1.))
    else:
        hh_path = txt_path.replace('.txt', f'_{n}p.png')
        hh, _ = ff.read_png(hh_path, alimit=pi2limit)
    capA, _ = gf.path_parts(hh_path)
    pf.plotAAB(hh, 'HHn', capA=capA, sxy=sxy, pause=.5)
    hhhp += [hh]
hha = hhhp[0]

lam_1ns = np.zeros(nh + 1)
lam1 = lam_ns[1]
for n in range(2, nh+1):
    lam_1ns[n] = (lam1 * lam_ns[n]) / (lam1 - lam_ns[n])
print(gf.prn_list('lam_1ns', lam_1ns, 1))

zz_1ns = [blank, blank]
for n in range(2, nh+1):
    zz1n = (np.mod(hhhp[1] - hhhp[n] + pi, pi2) - pi) * lam_1ns[n] / pi2
    zz_1ns += [zz1n]
    pf.plotAAB(zz1n, 'ZZ1n', capA=f'ZZ1{n}', sxy=sxy, pause=.5)

zz_12ns = [blank, blank, zz_1ns[2]]
for n in range(3, nh+1):
    zz_12n = df.stitch(zz_12ns[n-1], zz_1ns[n], lam_1ns[2], lam_1ns[n])
    zz_12ns += [zz_12n]
    pf.plotAAB(zz_12n, 'ZZ12n', capA=f'ZZ12{n}', sxy=sxy, pause=.5)

pf.plt.show()







