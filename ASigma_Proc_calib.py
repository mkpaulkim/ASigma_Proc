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
roi = (560, 950, 10, 10)

# txt_path = filedialog.askopenfilename(title='TXT file path', filetypes=[('txt files', '*.txt')])
txt_path = '/media/mkpaulkim/Ultra Touch/{{UT White}}/Dropbox/[[ PROJECTS.dbox ]]/project folders 2021/proj 2021-10 AlphaSigma/temp_data/aaa.txt'
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
    # pf.plotAAB(hh, capA=capA, roi=roi, sxy=sxy, pause=1)
    hhhp += [hh]
hha = hhhp[0]

lam_1ns = np.zeros(nh + 1)
lam1 = lam_ns[1]
for n in range(2, nh+1):
    lam_1ns[n] = (lam1 * lam_ns[n]) / (lam1 - lam_ns[n])
# lam_1ns[5] = 500
print(gf.prn_list('lam_1ns', lam_1ns, 1))
lam12 = lam_1ns[2]

zz_1ns = [blank, blank]
ep_1ns = [blank, blank]
for n in range(2, nh+1):
    ep1n = np.mod(hhhp[n] - hhhp[1] + pi, pi2) - pi
    zz1n = ep1n * lam_1ns[n] / pi2
    ep_1ns += [ep1n]
    zz_1ns += [zz1n]
    # pf.plotAAB(zz1n, capA=f'ZZ1{n}', roi=roi, sxy=sxy, pause=1)

zz_12ns = [blank, blank, zz_1ns[2]]
for n in range(3, nh+1):

    lam_1ns[n] = df.calib_lam1n(zz_12ns[n-1], ep_1ns[n], lam12, lam_1ns[n])
    zz_12n = df.stitch(zz_12ns[n-1], zz_1ns[n], lam_1ns[2], lam_1ns[n])
    zz_12ns += [zz_12n]
    pf.plotAAB(zz_12n, capA=f'ZZ12{n}: lam_1{n} = {lam_1ns[n]:.1f}', roi=roi, sxy=sxy, pause=1)

graphs = []
ix, iy, rx, ry = roi
print(f'< roi = {roi}')
for n in range(1, nh+1):
    graphs += [(hhhp[n][iy, :], (n-1, 0), f'HH{n}p', pi2limit)]
for n in range(2, nh+1):
    lam_1n = lam_1ns[n]
    graphs += [(zz_1ns[n][iy, :], (n-1, 1), f'ZZ_1{n}', (-lam_1n/2, lam_1n/2))]
for n in range(3, nh+1):
    graphs += [(zz_12ns[n][iy, :], (n-1, 2), f'ZZ_12{n}', (-lam12/2, lam12/2))]

pf.graph_many(graphs, col_row=(3, nh), sxy=(.25, .25), pause=1)

pf.plt.show()





