import os
import numpy as np
import pycubelib.tkinter_parts as tp
import pycubelib.files_functions as ff
import pycubelib.general_functions as gf
import pycubelib.plotting_functions as pf
import pycubelib.data_functions as df

# cwd = os.getcwd()
pi = np.pi
pi2 = pi * 2
pi2limit = (-pi, pi)
sxy = (1, 1)

fp = tp.tkwindow('AlphaSigmaProc', (20, 50, 1050, 400), tkbg='gray85')

btn_readtxt = tp.CmdButton(fp, (100, 50, 10), 'read TXT')
btn_readphs = tp.CmdButton(fp, (100, 100, 10), 'read HHp')
btn_lam1ns = tp.CmdButton(fp, (100, 150, 10), 'get lam_1ns')
btn_nextn = tp.CmdButton(fp, (100, 200, 10), 'next n')
btn_inv = tp.CmdButton(fp, (450, 200, 5), 'inv', 'gray90')
btn_makezz = tp.CmdButton(fp, (100, 250, 10), 'make ZZ_n')
btn_note = tp.CmdButton(fp, (100, 350, 10), 'notes')
btn_save = tp.CmdButton(fp, (250, 350, 10), 'save')
btn_detail = tp.CmdButton(fp, (450, 250, 5), 'detail', 'gray90')
btn_proc = tp.CmdButton(fp, (850, 200, 10), 'proc')
btn_graphall = tp.CmdButton(fp, (850, 250, 10), 'graph all')
btn_mayavi = tp.CmdButton(fp, (850, 300, 10), 'mayavi')
btn_adios = tp.CmdButton(fp, (850, 50, 10), 'adios', 'indian red')

ent_txtpath = tp.ParamEntry(fp, (220, 55, 35), '', '')
ent_phspath = tp.ParamEntry(fp, (220, 105, 35), '', '')
ent_n = tp.ParamEntry(fp, (250, 205, 5), 0, 'n')
ent_nw = tp.ParamEntry(fp, (350, 205, 5), 0, 'nw')
ent_lam1n = tp.ParamEntry(fp, (300, 250, 15), 0, 'lam1n')
ent_roi = tp.ParamEntry(fp, (600, 50, 25), '1200, 900, 10, 10', 'ixy rxy')
ent_z0roi = tp.ParamEntry(fp, (600, 100, 15), 800.0, 'Z0_roi')
ent_zh = tp.ParamEntry(fp, (600, 150, 15), 0.0, 'Zh')
ent_qxy = tp.ParamEntry(fp, (600, 200, 25), '-0.3, 0, -3e-5', 'qxys')
ent_mnf = tp.ParamEntry(fp, (600, 250, 15), '3, 1', 'mnf')

prog_n = tp.ProgressBar(fp, (100, 310, 410), '')

txt_path = ''
roi = (0, 0, 0, 0)
z0_roi = 0.
zh = 0.
qxys = (0., 0., 0.)
mnf = (0, 0)
txt_note = ''
notes = ''

nxydw = (0, 0, 0, 0)
wln = []
hhh = []
lam_1ns = []
blank = []
zz_ns = []
zz_1ns = []
zz_12ns = []
zz_proc = []
noise = []
proc_noise = 0.0


def program_loop():
    global txt_path, roi, z0_roi, zh, qxys, mnf

    txt_path = ent_txtpath.get_val(str)
    roi = ent_roi.get_list_val()
    z0_roi = ent_z0roi.get_val(float)
    zh = ent_zh.get_val(float)
    qxys = ent_qxy.get_list_val(float)
    mnf = ent_mnf.get_list_val()

    fp.after(tloop, program_loop)


def read_txt():
    global txt_path, nxydw, wln, blank, txt_note

    text, txt_path = ff.read_txt()
    nx = gf.find_param(text, 'nx')
    ny = gf.find_param(text, 'ny')
    dx = gf.find_param(text, 'dx', float)
    nw = gf.find_param(text, 'nw')
    wln = gf.find_param(text, 'wln', float)
    nxydw = (nx, ny, dx, nw)
    blank = np.zeros((ny, nx))

    ent_txtpath.set_entry(txt_path)
    ent_nw.set_entry(nw)
    ent_n.set_entry(0)

    txt_note = f'>> txt_path = {txt_path}'
    txt_note += f'\n{text}'

    if len(wln) < nw + 1:
        wln = [0.0] + wln
        txt_note += f'\n> !!! len(wln) < nw + 1 = {nw + 1}'
        txt_note += f'\n> new wln: {gf.prn_list("wl", wln)} \n'

    print(txt_note)
    print(f'> read TXT: done ...')


def read_phs():
    global hhh
    nx, ny, dx, nw = nxydw

    print()
    hhh = []
    for m in range(nw+1):
        if m == 0:
            hh_path = txt_path.replace('.txt', f'_aa.png')
            hh, _ = ff.read_png(hh_path, alimit=(0., 1.))
        else:
            hh_path = txt_path.replace('.txt', f'_{m}p.png')
            hh, _ = ff.read_png(hh_path, alimit=pi2limit)
        hhh += [hh]

        ent_n.set_entry(m)
        prog_n.setval(100 * m / nw)
        ent_phspath.set_entry(hh_path)
        capA, _ = gf.path_parts(hh_path)
        pf.plotAAB(hh, capA=capA, roi=roi, sxy=sxy)

        print(f'> hh_path = {hh_path}')

    print(f'> read HHp: done ...')


def get_lam1ns():
    global lam_1ns, zz_ns, zz_1ns, zz_12ns, noise
    nx, ny, dx, nw = nxydw

    lam_1ns = [0, 0]
    lam1 = wln[1]
    for m in range(2, nw+1):
        lam1n = (lam1 * wln[m]) / (lam1 - wln[m])
        lam_1ns += [lam1n]

    zz_ns = [blank] * (nw+1)
    zz_1ns = [blank] * (nw+1)
    zz_12ns = [blank] * (nw+1)
    noise = [0] * (nw + 1)

    ent_n.set_entry(0)
    prog_n.setval(0)
    capA = gf.path_parts(ent_txtpath.get_val(str))[0].replace('.txt', '_aa.png')
    pf.plotAAB(hhh[0], capA=capA, roi=roi, sxy=sxy, ulimit=(0, 1))

    print('\n>>> get_lam1ns: reset zz_ns, zz_1ns, zz_12ns ...')
    print(gf.prn_list('lam_1ns', lam_1ns, 1) + '\n')


def nextn():
    nx, ny, dx, nw = nxydw

    n = np.mod(ent_n.get_val() + 1, nw+1)
    ent_n.set_entry(n)
    prog_n.setval(100 * n / nw)
    ent_lam1n.set_entry(f'{lam_1ns[n]:.1f}')
    make_zz()


def make_zz():
    global zz_ns, zz_1ns, zz_12ns, zz_proc, noise

    m = ent_n.get_val()
    sign = 1 - btn_inv.is_on() * 2
    wl = wln[m]
    lam1n = lam_1ns[m] = ent_lam1n.get_val(float)
    lam12 = lam_1ns[2]
    hha = hhh[0].copy()
    hhp = hhh[m].copy()

    if m >= 1:
        hhp, hha = df.diffract(hhh[m], hha, wl, nxydw, zh)
        zz_ns[m] = sign * hhp * wl / pi2

        pf.plotAAB(zz_ns[m], figname='ZZn', capA=f'ZZ_{m}', capB=f'lam_{m} = {wl:.8f}',
                   roi=roi, sxy=sxy, ulimit=(-wl/2, wl/2))

    if m >= 2:
        zz1n = (np.mod((zz_ns[1] * pi2 / wln[1] - zz_ns[m] * pi2 / wln[m]) + pi, pi2) - pi) * lam1n / pi2
        z_roi, _ = df.roi_cyclic_measure(zz1n, roi, lam1n)
        zz_1ns[m] = np.mod(zz1n - z_roi + z0_roi + lam1n/2, lam1n) - lam1n/2

        pf.plotAAB(zz_1ns[m], figname='ZZ1n', capA=f'ZZ_1{m}', capB=f'lam_1{m} = {lam_1ns[m]:.1f}',
                   roi=roi, sxy=sxy, ulimit=(-lam1n/2, lam1n/2))

    if m == 2:
        zz_12ns[2] = zz_1ns[2].copy()

    if m >= 3:
        graphs, gxy = df.calib_lam1n(zz_12ns[m-1], zz_1ns[m], lam1n, roi)
        if btn_detail.is_on():
            pf.graph_many(graphs, 'calibrate', gxy)

        zz_12ns[m], graphs, gxy = df.stitch(zz_12ns[m-1], zz_1ns[m], lam12, lam1n, roi)
        if btn_detail.is_on():
            pf.graph_many(graphs, 'stitch', gxy, sxy=(1, .75))

    if m >= 2:
        zz_12ns[m] = df.cyclic_medfilter(zz_12ns[m], mnf, lam12)
        zz_proc = zz_12ns[m].copy()

        _, noise[m] = df.roi_measure(zz_12ns[m], roi)
        pf.plotAAB(zz_12ns[m], figname='ZZ12n', capA=f'ZZ_12{m}',
                   capB=f'lam12 = {lam_1ns[2]:.1f}; lam1{m} = {lam_1ns[m]:.1f}; noise = {noise[m]:.1f}',
                   roi=roi, sxy=sxy, ulimit=(-lam12/2, lam12/2))

    print(f'> make_zz: m = {m}; sign = {sign}; wl_{m} = {wl:.8f}; lam_1{m} = {lam1n:.1f}; ' 
          f'zh = {zh:.1f}; z0_roi = {z0_roi:.1f}; mnf = {mnf}')


def proc_zz():
    global zz_proc, proc_noise

    n = ent_n.get_val()
    lam12 = lam_1ns[2]
    zz_proc = df.zz_tilt(zz_12ns[n], nxydw, qxys, lam12)
    _, proc_noise = df.roi_measure(zz_proc, roi)
    pf.plotAAB(zz_proc, figname='ZZproc', capA=f'ZZ_proc', capB=f'lam12 = {lam_1ns[2]:.1f}; noise = {proc_noise:.1f}',
               roi=roi, sxy=sxy, ulimit=(-lam12/2, lam12/2))

    print(f'\n> proc_zz: qxys = {qxys}')


def mayavi():
    n = ent_n.get_val()
    cap = gf.path_parts(txt_path)[0] + f': ZZ_12{n}'
    lam12 = lam_1ns[2]
    pf.mayaviAA(zz_proc, caption=cap, ulimit=(-lam12/2, lam12/2))

    print(f'> mayavi: done ...')


def graph_all():
    nx, ny, dx, nw = nxydw

    lam12 = lam_1ns[2]
    graphs = []
    ix, iy, rx, ry = roi
    for m in range(1, nw + 1):
        wl = wln[m]
        graphs += [(zz_ns[m][iy, :], (m - 1, 0), f'ZZ{m}p: wl{m} = {wl:.8f}', (0, nx), (-wl/2, wl/2))]
    for m in range(2, nw + 1):
        lam_1n = lam_1ns[m]
        graphs += [(zz_1ns[m][iy, :], (m - 1, 1), f'ZZ_1{m}: lam1{m} = {lam_1n:.1f}', (0, nx), (-lam_1n/2, lam_1n/2))]
        graphs += [(zz_12ns[m][iy, :], (m - 1, 2),
                    f'ZZ_12{m}: lam_1{m} = {lam_1ns[m]: .1f}, noise = {noise[m]:.1f}', (0, nx), (-lam12/2, lam12/2))]
    graphs += [(zz_proc[iy, :], (0, 2), f'ZZ_proc', (0, nx), (-lam12/2, lam12/2))]

    pf.graph_many(graphs, col_row=(3, nw), sxy=(.75, .75), pause=1)

    print(f'> graph_all: done ...')


def print_notes():
    global notes

    notes = txt_note
    notes += f'> ================'
    notes += '\n' + gf.runstamp(__file__)
    notes += '\n' + gf.prn_list('lam_1ns', lam_1ns, 1)
    notes += '\n' + gf.prn_list('roi', roi, 0) + f'; z0_roi = {z0_roi:.1f}'
    notes += '\n' + gf.prn_list('qxys', qxys, 3) + f'; zh = {zh:.1f}'
    notes += '\n' + gf.prn_list('mnf', mnf, 0)
    notes += '\n' + gf.prn_list('noise', noise, 1)[:-2] + f' + [{proc_noise:.1f}]; '

    print()
    print(f'> notes ---------------------------------------------------------------------')
    print(notes)
    print(f'> ---------------------------------------------------------------------------')
    print()


def save_zzproc():
    from tkinter import filedialog

    proc_txt_file, _ = gf.path_parts(txt_path.replace('.txt', '_proc.txt'), 1)
    proc_txt_path = filedialog.asksaveasfilename(initialfile=proc_txt_file)
    png_path = proc_txt_path.replace('.txt', '.png')

    ff.write_txt(notes, proc_txt_path)
    ff.write_png(zz_proc, png_path, alimit=(-lam_1ns[2]/2, lam_1ns[2]/2))


def adios():
    # print(notes)
    fp.destroy()
    print(f'> adios amigos ...')
    quit()


btn_readtxt.command(read_txt)
btn_readphs.command(read_phs)
btn_lam1ns.command(get_lam1ns)
btn_nextn.command(nextn)
btn_makezz.command(make_zz)
btn_note.command(print_notes)
btn_save.command(save_zzproc)
btn_inv.command(btn_inv.switch)
btn_detail.command(btn_detail.switch)
btn_proc.command(proc_zz)
btn_graphall.command(graph_all)
btn_mayavi.command(mayavi)
btn_adios.command(adios)

btn_inv.on()

tloop = 10
fp.after(tloop, program_loop)
fp.mainloop()






