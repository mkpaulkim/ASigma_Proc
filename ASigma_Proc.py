import numpy as np
import pycubelib.tkinter_parts as tp
import pycubelib.files_functions as ff
import pycubelib.general_functions as gf
import pycubelib.plotting_functions as pf
import pycubelib.data_functions as df

pi = np.pi
pi2 = pi * 2
pi2limit = (-pi, pi)
sxy = (1, 1)

nxydw = (0, 0, 0, 0)
roi = (0, 0, 0, 0)
z0_roi = 0.0
wln = []
hhh = []
lam_1ns = []
zz_ns = []
zz_1ns = []
zz_12ns = []
zz_proc = []
noise = []

fp = tp.tkwindow('AlphaSigmaProc', (20, 50, 1050, 400), tkbg='gray90')

btn_readtxt = tp.CmdButton(fp, (100, 50, 10), 'read TXT', 'orange')
ent_txtpath = tp.ParamEntry(fp, (220, 55, 35), '', '')
btn_readphs = tp.CmdButton(fp, (100, 100, 10), 'read HHp', 'orange')
ent_phspath = tp.ParamEntry(fp, (220, 105, 35), '', '')
btn_lam1ns = tp.CmdButton(fp, (100, 150, 10), 'get lam_1ns', 'orange')
btn_nextn = tp.CmdButton(fp, (100, 200, 10), 'next n', 'orange')
ent_n = tp.ParamEntry(fp, (250, 205, 5), 0, 'n')
ent_nw = tp.ParamEntry(fp, (350, 205, 5), 0, 'nw')
btn_inv = tp.CmdButton(fp, (450, 200, 5), 'inv')
btn_makezz = tp.CmdButton(fp, (100, 250, 10), 'make ZZ_n', 'orange')
ent_lam1n = tp.ParamEntry(fp, (300, 250, 15), 0, 'lam_1n')
btn_detail = tp.CmdButton(fp, (450, 250, 5), 'detail')
prog_n = tp.ProgressBar(fp, (100, 310, 410), '')
ent_roi = tp.ParamEntry(fp, (600, 50, 25), '1200, 900, 10, 10', 'ixy rxy')
ent_z0 = tp.ParamEntry(fp, (600, 100, 15), 800.0, 'Z0_roi')
ent_zh = tp.ParamEntry(fp, (600, 150, 15), 0.0, 'Zh')
ent_qxy = tp.ParamEntry(fp, (600, 200, 25), '-0.3, 0, -3e-5', 'qxys')
ent_mnf = tp.ParamEntry(fp, (600, 250, 15), '3, 1', 'mnf')
btn_diffract = tp.CmdButton(fp, (850, 150, 10), 'diffract')
btn_tilt = tp.CmdButton(fp, (850, 200, 10), 'tilt', 'orange')
btn_graphall = tp.CmdButton(fp, (850, 250, 10), 'graph all', 'orange')
btn_mayavi = tp.CmdButton(fp, (850, 300, 10), 'mayavi', 'orange')

btn_adios = tp.CmdButton(fp, (850, 50, 10), 'adios', 'indian red')


def program_loop():
    global roi, z0_roi



    roi = tuple(ent_roi.get_list_val())
    z0_roi = ent_z0.get_val(float)

    fp.after(tloop, program_loop)


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

    if len(wln) < nw + 1:
        print(f'> !!! len(wln) = {len(wln)} < nw + 1 = {nw + 1}')
        wln = [0.0] + wln
        print(f'> new wln: {gf.prn_list("wln", wln)}')


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
    global lam_1ns, zz_ns, zz_1ns, zz_12ns, noise
    nx, ny, dx, nw = nxydw

    lam_1ns = [0, 0]
    lam1 = wln[1]
    for n in range(2, nw+1):
        lam1n = (lam1 * wln[n]) / (lam1 - wln[n])
        lam_1ns += [lam1n]

    zz_ns = [blank] * (nw+1)
    zz_1ns = [blank] * (nw+1)
    zz_12ns = [blank] * (nw+1)
    noise = [0] * (nw + 1)

    print(gf.prn_list('lam_1ns', lam_1ns, 1))
    ent_n.set_entry(0)
    prog_n.setval(0)
    capA = gf.path_parts(ent_txtpath.get_val(str))[0].replace('.txt', '_aa.png')
    pf.plotAAB(hhh[0], capA=capA, roi=roi, sxy=sxy)


def nextn():
    nx, ny, dx, nw = nxydw

    n = np.mod(ent_n.get_val() + 1, nw+1)
    ent_n.set_entry(n)
    prog_n.setval(100 * n / nw)
    ent_lam1n.set_entry(f'{lam_1ns[n]:.1f}')
    make_zz()


def make_zz():
    global zz_ns, zz_1ns, zz_12ns, zz_proc, noise

    n = ent_n.get_val()
    sign = 1 - btn_inv.is_on() * 2
    lam1n = lam_1ns[n] = ent_lam1n.get_val(float)
    lam12 = lam_1ns[2]
    hha = hhh[0].copy()
    hhp = hhh[n].copy()

    if n >= 1:
        if btn_diffract.is_on():
            zh = tuple(ent_zh.get_list_val(float))
            hhp, hha = df.diffract(hhh[n], hha, wln[n], nxydw, zh)
        zz_ns[n] = sign * hhp * wln[n] / pi2

        pf.plotAAB(zz_ns[n], figname='ZZn', capA=f'ZZ_{n}', capB=f'lam_{n} = {wln[n]:.8f}', roi=roi, sxy=sxy)

    if n >= 2:
        zz1n = (np.mod((zz_ns[1] * pi2 / wln[1] - zz_ns[n] * pi2 / wln[n]) + pi, pi2) - pi) * lam1n / pi2
        z_roi, _ = df.roi_cyclic_measure(zz1n, roi, lam1n)
        zz_1ns[n] = np.mod(zz1n - z_roi + z0_roi + lam1n/2, lam1n) - lam1n/2

        pf.plotAAB(zz_1ns[n], figname='ZZ1n', capA=f'ZZ_1{n}',
                   capB=f'lam_1{n} = {lam_1ns[n]:.1f}', roi=roi, sxy=sxy)

    if n == 2:
        zz_12ns[2] = zz_1ns[2].copy()

    if n >= 3:
        graphs, gxy = df.calib_lam1n(zz_12ns[n-1], zz_1ns[n], lam1n, roi)
        if btn_detail.is_on():
            pf.graph_many(graphs, 'calibrate', gxy)

        zz_12ns[n], graphs, gxy = df.stitch(zz_12ns[n-1], zz_1ns[n], lam12, lam1n, roi)
        if btn_detail.is_on():
            pf.graph_many(graphs, 'stitch', gxy, sxy=(1, .75))

    if n >= 2:
        if btn_tilt.is_on():
            qxys = tuple(ent_qxy.get_list_val(float))
            zz_12ns[n] = df.zz_tilt(zz_12ns[n], nxydw, qxys, lam_1ns[2])
            btn_tilt.off()

        mnf = tuple(ent_mnf.get_list_val())
        zz_12ns[n] = df.cyclic_medfilter(zz_12ns[n], mnf, lam12)
        zz_proc = zz_12ns[n].copy()

        _, noise[n] = df.roi_measure(zz_12ns[n], roi)
        pf.plotAAB(zz_12ns[n], figname='ZZ12n', capA=f'ZZ_12{n}', roi=roi, sxy=sxy,
                   capB=f'lam12 = {lam_1ns[2]:.1f}; lam1{n} = {lam_1ns[n]:.1f}; noise = {noise[n]:.1f}')


def make_zz_proc():
    global zz_proc

    n = ent_n.get_val()
    qxys = tuple(ent_qxy.get_list_val(float))
    lam12 = lam_1ns[2]
    zz_proc = df.zz_tilt(zz_12ns[n], nxydw, qxys, lam12)
    _, noise[n] = df.roi_measure(zz_proc, roi)
    pf.plotAAB(zz_proc, figname='ZZ12n', capA=f'ZZ_proc', roi=roi, sxy=sxy,
               capB=f'lam12 = {lam_1ns[2]:.1f}; noise = {noise[n]:.1f}')


def mayavi():
    n = ent_n.get_val()
    cap = gf.path_parts(ent_txtpath.get_val(str))[0] + f': ZZ_12{n}'
    pf.mayaviAA(zz_proc, caption=cap)


def graph_all():
    nx, ny, dx, nw = nxydw

    lam12 = lam_1ns[2]
    graphs = []
    ix, iy, rx, ry = roi
    for n in range(1, nw + 1):
        wl = wln[n]
        graphs += [(zz_ns[n][iy, :], (n - 1, 0), f'ZZ{n}p: wl{n} = {wl:.8f}', (), (-wl/2, wl/2))]
    for n in range(2, nw + 1):
        lam_1n = lam_1ns[n]
        graphs += [(zz_1ns[n][iy, :], (n - 1, 1), f'ZZ_1{n}: lam1{n} = {lam_1n:.1f}', (), (-lam_1n/2, lam_1n/2))]
        graphs += [(zz_12ns[n][iy, :], (n - 1, 2),
                    f'ZZ_12{n}: lam_1{n} = {lam_1ns[n]: .1f}, noise = {noise[n]:.1f}', (), (-lam12/2, lam12/2))]
    graphs += [(zz_proc[iy, :], (0, 2), f'ZZ_proc', (), (-lam12/2, lam12/2))]

    pf.graph_many(graphs, col_row=(3, nw), sxy=(.75, .75), pause=1)


def adios():
    fp.destroy()
    print(f'> adios amigos ...')
    quit()


btn_readtxt.command(read_txt)
btn_readphs.command(read_phs)
btn_lam1ns.command(get_lam1ns)
btn_nextn.command(nextn)
btn_makezz.command(make_zz)
btn_inv.command(btn_inv.switch)
btn_detail.command(btn_detail.switch)
btn_diffract.command(btn_diffract.switch)
btn_tilt.command(make_zz_proc)
btn_graphall.command(graph_all)
btn_mayavi.command(mayavi)
btn_adios.command(adios)

btn_inv.on()

tloop = 10
fp.after(tloop, program_loop)
fp.mainloop()






