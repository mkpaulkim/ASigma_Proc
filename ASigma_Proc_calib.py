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
roi = (1000, 900, 10, 10)
z0_roi = 0.

# txt_path = filedialog.askopenfilename(title='TXT file path', filetypes=[('txt files', '*.txt')])
txt_path = '/media/mkpaulkim/Ultra Touch/{{UT White}}/Dropbox/[[ PROJECTS.dbox ]]' \
           '/project folders 2021/proj 2021-10 AlphaSigma/temp_data/eee7.txt'
notes, txt_path = ff.read_txt(txt_path)
note = notes[notes.find('%%%'):]
nx = gf.find_param(note, 'nx', int)
ny = gf.find_param(note, 'ny', int)
dx = gf.find_param(note, 'dx', float)
nw = gf.find_param(note, 'nw', int)
wln = gf.find_param(note, 'wln', float)
if wln[0] != 0.:
    wln = [0.] + wln
blank = np.zeros((ny, nx))

print(f'> txt_path = {txt_path}')
print(f'> notes: \n{notes}')

''' read hhh '''
hhhp = []
for n in range(nw + 1):
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

''' get lam_1ns '''
lam_1ns = np.zeros(nw + 1)
lam1 = wln[1]
for n in range(2, nw + 1):
    lam_1ns[n] = (lam1 * wln[n]) / (lam1 - wln[n])
# lam_1ns[5] = 500
lam12 = lam_1ns[2]

print(gf.prn_list('wln', wln, 8))
print(gf.prn_list('lam_1ns', lam_1ns, 1))

''' make zz_1ns '''
lam_1ns0 = lam_1ns.copy()
ep_1ns = [blank, blank]
zz_1ns = [blank, blank]
ave_1ns = [0, 0]
std_1ns = [0, 0]
for n in range(2, nw + 1):
    lam1n = lam_1ns[n]
    ep1n = np.mod(hhhp[n] - hhhp[1] + pi, pi2) - pi
    zp1n = ep1n * lam1n / pi2
    z_roi, _ = gf.roi_cyclic_measure(zp1n, roi, lam1n)
    zz1n = np.mod(zp1n - z_roi + z0_roi + lam1n/2, lam1n) - lam1n/2

    ave1n, std1n = gf.roi_measure(zz1n, roi)
    # _, z1_roi = gf.roi_measure(zz1n, roi)
    # print(f'< z_roi = {z_roi:.1f}, z0_roi = {z0_roi:.1f}, z1_roi = {z1_roi:.1f}')

    ep_1ns += [ep1n]
    zz_1ns += [zz1n]
    ave_1ns += [ave1n]
    std_1ns += [std1n]

    pf.plotAAB(zz1n, capA=f'ZZ1{n}', capB=f'ave = {ave1n:.1f}; std = {std1n:.1f}', roi=roi, sxy=sxy, pause=1)

print(gf.prn_list('ave_1ns', ave_1ns, 1))
print(gf.prn_list('std_1ns', std_1ns, 1))
pf.plt.close()

''' make zz_12ns '''
zz_12ns = [blank, blank, zz_1ns[2]]
ave_12ns = [0, 0, ave_1ns[2]]
std_12ns = [0, 0, std_1ns[2]]
for n in range(3, nw + 1):
    lam_1ns[n] = df.calib_lam1n(zz_12ns[n-1], zz_1ns[n], lam12, lam_1ns[n], roi)
    zz12n = df.stitch(zz_12ns[n-1], zz_1ns[n], lam12, lam_1ns[n])

    ave12n, std12n = gf.roi_measure(zz12n, roi)

    zz_12ns += [zz12n]
    ave_12ns += [ave12n]
    std_12ns += [std12n]

    pf.plotAAB(zz12n, capA=f'ZZ12{n}: lam_1{n} = {lam_1ns[n]:.1f}', capB=f'ave = {ave12n:.1f}; std = {std12n:.1f}', roi=roi, sxy=sxy, pause=1)

print(gf.prn_list('lam_1ns_old', lam_1ns0, 1))
print(gf.prn_list('lam_1ns_new', lam_1ns, 1))
print(gf.prn_list('ave_12ns', ave_12ns, 1))
print(gf.prn_list('std_12ns', std_12ns, 1))

''' graph all '''
graphs = []
ix, iy, rx, ry = roi
# print(f'< roi = {roi}')
for n in range(1, nw + 1):
    graphs += [(hhhp[n][iy, :], (n-1, 0), f'HH{n}p', pi2limit)]
for n in range(2, nw + 1):
    lam_1n = lam_1ns[n]
    graphs += [(zz_1ns[n][iy, :], (n-1, 1), f'ZZ_1{n}', (-lam_1n/2, lam_1n/2))]
for n in range(3, nw + 1):
    graphs += [(zz_12ns[n][iy, :], (n-1, 2), f'ZZ_12{n}', (-lam12/2, lam12/2))]

pf.graph_many(graphs, col_row=(3, nw), sxy=(.25, .25), pause=1)

pf.mayaviAA(zz_12ns[-1])

pf.plt.show()





