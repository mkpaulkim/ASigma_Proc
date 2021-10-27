import numpy as np
from tkinter import filedialog
import pycubelib.front_panel as fp
import pycubelib.files_functions as ff
import pycubelib.general_functions as gf
import pycubelib.plotting_functions as pf
import pycubelib.data_functions as df

pi = np.pi
pi2 = pi * 2
pi2limit = (-pi, pi)

sxy = (.35, .35)
roi = (1000, 900, 10, 10)
z0_roi = 0

nxydw = (0, 0, 0, 0)
wln = []
hhh = []
lam_1ns = []
zz_ns = []
zz_1ns = []
zz_12ns = []

tkw = fp.tkwindow('AlphaSigmaProc', (20, 20, 1000, 350), tkbg='gray90')

btn_readtxt = fp.CmdButton(tkw, (100, 30, 10), 'read TXT', 'orange')
ent_txtpath = fp.ParamEntry(tkw, (220, 35, 35), '', '')
btn_readphs = fp.CmdButton(tkw, (100, 60, 10), 'read HHp', 'orange')
ent_phspath = fp.ParamEntry(tkw, (220, 65, 35), '', '')
btn_lam1ns = fp.CmdButton(tkw, (100, 90, 10), 'get lam_1ns', 'orange')
ent_n = fp.ParamEntry(tkw, (250, 155, 5), 0, 'n')
ent_nw = fp.ParamEntry(tkw, (350, 155, 5), 1, 'nw', 'r')
btn_zz1n = fp.CmdButton(tkw, (100, 150, 10), 'ZZ_1n', 'orange')

prog_readphs = fp.ProgressBar(tkw, (100, 210, 410), '')


def program_loop():

    # prog_readphs.setval(100 * ent_n.get_val() / ent_nw.get_val())

    tkw.after(tloop, program_loop)


def read_txt():
    global nxydw, wln, blank

    text, txt_path = ff.read_txt()
    nx = gf.find_param(text, 'nx')
    ny = gf.find_param(text, 'ny')
    dx = gf.find_param(text, 'dx', float)
    nw = gf.find_param(text, 'nw')
    wln = gf.find_param(text, 'wln', float)
    nxydw = (nx, ny, dx, nw)
    blank = np.zeros((ny, nx))

    ent_txtpath.set_entry(txt_path)
    print(f'>> txt_path = {txt_path} \n{text}')
    ent_nw.set_entry(nw)
    ent_n.set_entry(0)


def read_phs():
    global hhh
    nx, ny, dx, nw = nxydw

    txt_path = ent_txtpath.get_val(str)
    hhh = []
    for n in range(nw+1):
        if n == 0:
            hh_path = txt_path.replace('.txt', f'_aa.png')
            hh, _ = ff.read_png(hh_path, alimit=(0., 1.))
        else:
            hh_path = txt_path.replace('.txt', f'_{n}p.png')
            hh, _ = ff.read_png(hh_path, alimit=pi2limit)
        hhh += [hh]

        ent_n.set_entry(n)
        prog_readphs.setval(100 * (n) / nw)
        ent_phspath.set_entry(hh_path)
        capA, _ = gf.path_parts(hh_path)
        pf.plotAAB(hh, capA=capA, roi=roi, sxy=sxy, pause=1)


def get_lam1ns():
    global lam_1ns
    nx, ny, dx, nw = nxydw

    lam_1ns = [0, 0]
    lam1 = wln[1]
    for n in range(2, nw+1):
        lam1n = (lam1 * wln[n]) / (lam1 - wln[n])
        lam_1ns += [lam1n]

    print(gf.prn_list('lam_1ns', lam_1ns))
    ent_n.set_entry(0)
    prog_readphs.setval(0)
    pf.plt.close('plotAAB')


def make_zz1ns():
    global zz_ns, zz_1ns, zz_12ns
    nx, ny, dx, nw = nxydw

    n = ent_n.get_val()
    prog_readphs.setval(100 * n / nw)

    if n == 0:
        zz_ns = [blank] * nw
        zz_1ns = [blank] * nw
        zz_12ns = [blank] * nw

    if n >= 1:
        zz_ns[n] = hhh[n] * wln[n] / pi2

        pf.plotAAB(zz_ns[n], figname='ZZn', capA=f'ZZ_{n}', roi=roi, sxy=sxy, pause=1)

    if n >= 2:
        lam1n = lam_1ns[n]
        ep1n = np.mod(hhh[n] - hhh[1] + pi, pi2) - pi
        zz1n_ = ep1n * lam1n / pi2

        z_roi, _ = gf.roi_cyclic_measure(zz1n_, roi, lam1n)
        zz_1ns[n] = np.mod(zz1n_ - z_roi + z0_roi + lam1n/2, lam1n) - lam1n/2

        pf.plotAAB(zz_1ns[n], figname='ZZ1n', capA=f'ZZ_1{n}', roi=roi, sxy=sxy, pause=1)

    if n >= 3:
        lam12 = lam_1ns[2]
        # lam_1ns[n] = df.calib_lam1n(zz_12ns[n-1], zz_1ns[n], lam12, lam_1ns[n], roi)

        zz_12ns[n] = df.stitch(zz_12ns[n-1], zz_1ns[n], lam12, lam_1ns[n])

        pf.plotAAB(zz_12ns[n], figname='ZZ12n', capA=f'ZZ_12{n}', roi=roi, sxy=sxy, pause=1)

    n += 1
    if n > nw:
        n = 0
    ent_n.set_entry(n)


btn_readtxt.command(read_txt)
btn_readphs.command(read_phs)
btn_lam1ns.command(get_lam1ns)
btn_zz1n.command(make_zz1ns)

tloop = 10
tkw.after(tloop, program_loop)
tkw.mainloop()







#
# txt_path = filedialog.askopenfilename(title='TXT file path', filetypes=[('txt files', '*.txt')])
# # txt_path = '/media/mkpaulkim/Ultra Touch/{{UT White}}/Dropbox/[[ PROJECTS.dbox ]]' \
# #            '/project folders 2021/proj 2021-10 AlphaSigma/temp_data/eee7.txt'
# notes, txt_path = ff.read_txt(txt_path)
# note = notes[notes.find('%%%'):]
# nx = gf.find_param(note, 'nx', int)
# ny = gf.find_param(note, 'ny', int)
# dx = gf.find_param(note, 'dx', float)
# nw = gf.find_param(note, 'nw', int)
# wln = gf.find_param(note, 'wln', float)
# if wln[0] != 0.:
#     wln = [0.] + wln
# blank = np.zeros((ny, nx))
#
# print(f'> txt_path = {txt_path}')
# print(f'> notes: \n{notes}')
#
# ''' read hhh '''
# hhhp = []
# for n in range(nw + 1):
#     if n == 0:
#         hh_path = txt_path.replace('.txt', f'_aa.png')
#         hh, _ = ff.read_png(hh_path, alimit=(0., 1.))
#     else:
#         hh_path = txt_path.replace('.txt', f'_{n}p.png')
#         hh, _ = ff.read_png(hh_path, alimit=pi2limit)
#     capA, _ = gf.path_parts(hh_path)
#     # pf.plotAAB(hh, capA=capA, roi=roi, sxy=sxy, pause=1)
#     hhhp += [hh]
# hha = hhhp[0]
#
# ''' get lam_1ns '''
# lam_1ns = np.zeros(nw + 1)
# lam1 = wln[1]
# for n in range(2, nw + 1):
#     lam_1ns[n] = (lam1 * wln[n]) / (lam1 - wln[n])
# # lam_1ns[5] = 500
# lam12 = lam_1ns[2]
#
# print(gf.prn_list('wln', wln, 8))
# print(gf.prn_list('lam_1ns', lam_1ns, 1))
#
# ''' make zz_1ns '''
# lam_1ns0 = lam_1ns.copy()
# ep_1ns = [blank, blank]
# zz_1ns = [blank, blank]
# ave_1ns = [0, 0]
# std_1ns = [0, 0]
# for n in range(2, nw + 1):
#     lam1n = lam_1ns[n]
#     # ep1n = np.mod(hhhp[n] - hhhp[1] + pi, pi2) - pi
#     ep1n = np.mod(hhhp[1] - hhhp[n] + pi, pi2) - pi
#     zp1n = ep1n * lam1n / pi2
#     z_roi, _ = gf.roi_cyclic_measure(zp1n, roi, lam1n)
#     zz1n = np.mod(zp1n - z_roi + z0_roi + lam1n/2, lam1n) - lam1n/2
#
#     ave1n, std1n = gf.roi_measure(zz1n, roi)
#     # _, z1_roi = gf.roi_measure(zz1n, roi)
#     # print(f'< z_roi = {z_roi:.1f}, z0_roi = {z0_roi:.1f}, z1_roi = {z1_roi:.1f}')
#
#     ep_1ns += [ep1n]
#     zz_1ns += [zz1n]
#     ave_1ns += [ave1n]
#     std_1ns += [std1n]
#
#     pf.plotAAB(zz1n, capA=f'ZZ1{n}', capB=f'ave = {ave1n:.1f}; std = {std1n:.1f}', roi=roi, sxy=sxy, pause=1)
#
# print(gf.prn_list('ave_1ns', ave_1ns, 1))
# print(gf.prn_list('std_1ns', std_1ns, 1))
# pf.plt.close()
#
# ''' make zz_12ns '''
# zz_12ns = [blank, blank, zz_1ns[2]]
# ave_12ns = [0, 0, ave_1ns[2]]
# std_12ns = [0, 0, std_1ns[2]]
# for n in range(3, nw + 1):
#     lam_1ns[n] = df.calib_lam1n(zz_12ns[n-1], zz_1ns[n], lam12, lam_1ns[n], roi)
#     zz12n = df.stitch(zz_12ns[n-1], zz_1ns[n], lam12, lam_1ns[n])
#
#     ave12n, std12n = gf.roi_measure(zz12n, roi)
#
#     zz_12ns += [zz12n]
#     ave_12ns += [ave12n]
#     std_12ns += [std12n]
#
#     pf.plotAAB(zz12n, capA=f'ZZ12{n}: lam_1{n} = {lam_1ns[n]:.1f}', capB=f'ave = {ave12n:.1f}; std = {std12n:.1f}', roi=roi, sxy=sxy, pause=1)
#
# print(gf.prn_list('lam_1ns_old', lam_1ns0, 1))
# print(gf.prn_list('lam_1ns_new', lam_1ns, 1))
# print(gf.prn_list('ave_12ns', ave_12ns, 1))
# print(gf.prn_list('std_12ns', std_12ns, 1))
#
# ''' graph all '''
# graphs = []
# ix, iy, rx, ry = roi
# print(f'< roi = {roi}')
# for n in range(1, nw + 1):
#     graphs += [(hhhp[n][iy, :], (n-1, 0), f'HH{n}p', pi2limit)]
# for n in range(2, nw + 1):
#     lam_1n = lam_1ns[n]
#     graphs += [(zz_1ns[n][iy, :], (n-1, 1), f'ZZ_1{n}', (-lam_1n/2, lam_1n/2))]
# for n in range(3, nw + 1):
#     graphs += [(zz_12ns[n][iy, :], (n-1, 2), f'ZZ_12{n}', (-lam12/2, lam12/2))]
#
# pf.graph_many(graphs, col_row=(3, nw), sxy=(.25, .25), pause=1)
#
# pf.mayaviAA(zz_12ns[-1])
#
# pf.plt.show()
#
#
#


