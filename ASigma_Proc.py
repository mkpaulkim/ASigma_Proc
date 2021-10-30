import numpy as np
import pycubelib.front_panel as fp
import pycubelib.files_functions as ff
import pycubelib.general_functions as gf
import pycubelib.plotting_functions as pf
import pycubelib.data_functions as df

pi = np.pi
pi2 = pi * 2
pi2limit = (-pi, pi)
sxy = (1, 1)

nxydw = (0, 0, 0, 0)
wln = []
hhh = []
lam_1ns = []
zz_ns = []
zz_1ns = []
zz_12ns = []

tkw = fp.tkwindow('AlphaSigmaProc', (20, 50, 1000, 350), tkbg='gray90')

btn_readtxt = fp.CmdButton(tkw, (100, 50, 10), 'read TXT', 'orange')
ent_txtpath = fp.ParamEntry(tkw, (220, 55, 35), '', '')
btn_readphs = fp.CmdButton(tkw, (100, 100, 10), 'read HHp', 'orange')
ent_phspath = fp.ParamEntry(tkw, (220, 105, 35), '', '')
btn_lam1ns = fp.CmdButton(tkw, (100, 150, 10), 'get lam_1ns', 'orange')
ent_n = fp.ParamEntry(tkw, (250, 205, 5), 0, 'n')
ent_nw = fp.ParamEntry(tkw, (350, 205, 5), 0, 'nw', 'r')
btn_inv = fp.CmdButton(tkw, (450, 200, 5), 'inv')
btn_nextn = fp.CmdButton(tkw, (100, 200, 10), 'next n', 'orange')
btn_makezz = fp.CmdButton(tkw, (100, 250, 10), 'make ZZ_n', 'orange')
ent_lam1n = fp.ParamEntry(tkw, (300, 250, 15), 0, 'lam_1n')
btn_detail = fp.CmdButton(tkw, (450, 250, 5), 'detail')
prog_n = fp.ProgressBar(tkw, (100, 300, 410), '')
ent_roi = fp.ParamEntry(tkw, (700, 50, 25), '1000, 500, 10, 10', 'roi')
ent_z0 = fp.ParamEntry(tkw, (700, 100, 15), 0.0, 'Z0_roi')
ent_mnf = fp.ParamEntry(tkw, (700, 150, 15), '3, 1', 'mnf')

btn_adios = fp.CmdButton(tkw, (700, 300, 10), 'adios', 'indian red')


def program_loop():
    global roi, z0_roi
    nx, ny, dx, nw = nxydw

    roi = tuple(ent_roi.get_list_val())
    z0_roi = ent_z0.get_val(float)

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
        prog_n.setval(100 * (n) / nw)
        ent_phspath.set_entry(hh_path)
        capA, _ = gf.path_parts(hh_path)
        pf.plotAAB(hh, capA=capA, roi=roi, sxy=sxy)


def get_lam1ns():
    global lam_1ns, zz_ns, zz_1ns, zz_12ns
    nx, ny, dx, nw = nxydw

    lam_1ns = [0, 0]
    lam1 = wln[1]
    for n in range(2, nw+1):
        lam1n = (lam1 * wln[n]) / (lam1 - wln[n])
        lam_1ns += [lam1n]

    zz_ns = [blank] * (nw+1)
    zz_1ns = [blank] * (nw+1)
    zz_12ns = [blank] * (nw+1)

    print(gf.prn_list('lam_1ns', lam_1ns, 1))
    ent_n.set_entry(0)
    prog_n.setval(0)
    pf.plt.close('plotAAB')


def nextn():
    nx, ny, dx, nw = nxydw
    n = np.mod(ent_n.get_val() + 1, nw+1)
    ent_n.set_entry(n)
    prog_n.setval(100 * n / nw)
    ent_lam1n.set_entry(f'{lam_1ns[n]:.1f}')
    make_zz()


def make_zz():
    global zz_ns, zz_1ns, zz_12ns

    n = ent_n.get_val()
    sign = 1 - btn_inv.is_on() * 2
    lam1n = lam_1ns[n] = ent_lam1n.get_val(float)
    lam12 = lam_1ns[2]

    if n >= 1:
        zz_ns[n] = sign * hhh[n] * wln[n] / pi2

        pf.plotAAB(zz_ns[n], figname='ZZn', capA=f'ZZ_{n}', capB=f'lam_{n} = {wln[n]:.8f}', roi=roi, sxy=sxy)

    if n >= 2:
        ep1n = np.mod(sign * (hhh[1] - hhh[n]) + pi, pi2) - pi
        zz1n_ = ep1n * lam1n / pi2

        z_roi, _ = gf.roi_cyclic_measure(zz1n_, roi, lam1n)
        zz_1ns[n] = np.mod(zz1n_ - z_roi + z0_roi + lam1n/2, lam1n) - lam1n/2

        pf.plotAAB(zz_1ns[n], figname='ZZ1n', capA=f'ZZ_1{n}',
                   capB=f'lam_1{n} = {lam_1ns[n]:.1f}', roi=roi, sxy=sxy)

    if n == 2:
        zz_12ns[2] = zz_1ns[2].copy()

    if n >= 3:
        graphs, gxy = df.calib_lam1n(zz_12ns[n-1], zz_1ns[n], lam1n, roi)
        if btn_detail.is_on():
            pf.graph_many(graphs, 'calibrate', gxy)

        zz_12ns[n], graphs, gxy = df.stitch(zz_12ns[n-1], zz_1ns[n], lam12, lam_1ns[n], roi)
        if btn_detail.is_on():
            pf.graph_many(graphs, 'stitch', gxy, sxy=(1, .75))

    if n >= 2:
        mf, nf = tuple(ent_mnf.get_list_val())
        zz_proc = zz_12ns[n].copy()
        for i in range(nf):
            pass
            # zz_proc = np.med
        _, noise = gf.roi_measure(zz_12ns[n], roi)

        pf.plotAAB(zz_12ns[n], figname='ZZ12n', capA=f'ZZ_12{n}', roi=roi, sxy=sxy,
                   capB=f'lam12 = {lam_1ns[2]:.1f}; lam1{n} = {lam_1ns[n]:.1f}; noise = {noise:.1f}')


def adios():
    tkw.destroy()
    print(f'> adios amigos ...')
    quit()


btn_readtxt.command(read_txt)
btn_readphs.command(read_phs)
btn_lam1ns.command(get_lam1ns)
btn_nextn.command(nextn)
btn_makezz.command(make_zz)
btn_inv.command(btn_inv.switch)
btn_detail.command(btn_detail.switch)
btn_adios.command(adios)

tloop = 10
tkw.after(tloop, program_loop)
tkw.mainloop()






