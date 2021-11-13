# import os
import numpy as np
import pycubelib.tkinter_parts as tp
import pycubelib.files_functions as ff
import pycubelib.general_functions as gf
import pycubelib.plotting_functions as pf
import pycubelib.data_functions as df

pi = np.pi
pi2 = pi * 2
pi2limit = (-pi, pi)
sxy = (.85, 1)

fp = tp.tkwindow('AlphaSigmaProc', (20, 50, 1500, 400), tkbg='gray85')

btn_readtxt = tp.CmdButton(fp, (100, 50, 10), 'read TXT')
btn_readphs = tp.CmdButton(fp, (100, 100, 10), 'read HHp')
btn_lam1ns = tp.CmdButton(fp, (100, 150, 10), 'get lam_1ns')
btn_nextn = tp.CmdButton(fp, (100, 200, 10), 'next n')
btn_inv = tp.CmdButton(fp, (450, 200, 5), 'inv', 'gray90')
btn_makezz = tp.CmdButton(fp, (100, 250, 10), 'make ZZ_n')
btn_detail = tp.CmdButton(fp, (450, 250, 5), 'detail', 'gray90')
btn_graphall = tp.CmdButton(fp, (850, 250, 10), 'graph all')
btn_note = tp.CmdButton(fp, (850, 150, 10), 'notes')
btn_save = tp.CmdButton(fp, (850, 200, 10), 'save')
btn_plot = tp.CmdButton(fp, (850, 300, 10), 'plot')
btn_adios = tp.CmdButton(fp, (850, 50, 10), 'adios', 'indian red')

ent_txtpath = tp.ParamEntry(fp, (220, 55, 35), '', '')
ent_phspath = tp.ParamEntry(fp, (220, 105, 35), '', '')
ent_n = tp.ParamEntry(fp, (250, 205, 5), 0, 'n')
ent_nw = tp.ParamEntry(fp, (350, 205, 5), 0, 'nw')
ent_lam1n = tp.ParamEntry(fp, (300, 250, 15), 0, 'lam1n')
ent_roi = tp.ParamEntry(fp, (600, 50, 25), '1000, 500, 10, 10', 'ixy rxy')
ent_z0roi = tp.ParamEntry(fp, (600, 100, 15), 0.0, 'Z0_roi')
ent_zh = tp.ParamEntry(fp, (600, 150, 15), 0.0, 'Zh')
ent_mnf = tp.ParamEntry(fp, (600, 250, 15), '3, 1', 'mnf')
ent_np = tp.ParamEntry(fp, (800, 305, 5), 0, '')

prog_n = tp.ProgressBar(fp, (100, 310, 410), '')

btn_proc = tp.CmdButton(fp, (1100, 50, 10), 'proc')
ent_qxy = tp.ParamEntry(fp, (1100, 100, 25), '-0.300, 0.000, -0.020', 'qxys') # qx, qy in rad; ss curvature in 1/mm
ent_mnfp = tp.ParamEntry(fp, (1100, 150, 15), '3, 1', 'mnfp')
btn_dnoise = tp.CmdButton(fp, (1100, 200, 10), 'denoise')
ent_gamma = tp.ParamEntry(fp, (1100, 250, 15), 1.0, 'gamma')
btn_mayavi = tp.CmdButton(fp, (1100, 300, 10), 'mayavi')
btn_clear = tp.CmdButton(fp, (1300, 300, 10), 'clear')

txt_path = ''
roi = (0, 0, 0, 0)
# z0_roi = 0.
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
    global txt_path, roi, z0_roi, zh, mnf, qxys, mnfp

    txt_path = ent_txtpath.get_val(str)
    roi = ent_roi.get_list_val()
    z0_roi = ent_z0roi.get_val(float)
    zh = ent_zh.get_val(float)
    mnf = ent_mnf.get_list_val()
    qxys = ent_qxy.get_list_val(float)
    mnfp = ent_mnfp.get_list_val()

    fp.after(tloop, program_loop)


def read_txt():
    global txt_path, nxydw, wln, blank, txt_note

    text, txt_path = ff.read_txt()
    t = text.find('%%%')
    param_txt = text[t:]
    param_txt.replace('([', '[')            # 2018-01 Saturn data
    nx = gf.find_param(param_txt, 'nx')
    ny = gf.find_param(param_txt, 'ny')
    dx = gf.find_param(param_txt, 'dx', float)
    nw = gf.find_param(param_txt, 'nw')
    wln = gf.find_param(param_txt, 'wln', float)
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
    # print(f'< len(wln) = {len(wln)}: {wln}')
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
        pf.plotAAB(hh, figname='HHn', capA=capA, roi=roi, sxy=sxy)

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
    # capA = gf.path_parts(ent_txtpath.get_val(str))[0].replace('.txt', '_aa.png')
    capA = gf.path_parts(txt_path)[0].replace('.txt', '_aa.png')
    pf.plotAAB(hhh[0], figname='HHn', capA=capA, roi=roi, sxy=sxy, ulimit=(0, 1))

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
    nw = ent_nw.get_val()
    # z0_roi = ent_z0roi.get_val(float)
    sign = 1 - btn_inv.is_on() * 2
    wl = wln[m]
    lam1n = lam_1ns[m] = ent_lam1n.get_val(float)
    lam12 = lam_1ns[2]
    hha = hhh[0].copy()

    if 1 <= m <= nw:
        hhp, hha = df.diffract(hhh[m], hha, wl, nxydw, zh)
        zz_ns[m] = sign * hhp * wl / pi2

    if 2 <= m <= nw:
        zz1n = (np.mod((zz_ns[1] * pi2 / wln[1] - zz_ns[m] * pi2 / wln[m]) + pi, pi2) - pi) * lam1n / pi2
        z_roi, _ = df.roi_cyclic_measure(zz1n, roi, lam1n)
        zz_1ns[m] = np.mod(zz1n - z_roi + z0_roi + lam1n/2, lam1n) - lam1n/2

    if m == 2:
        zz_12ns[2] = zz_1ns[2].copy()

    if 3 <= m <= nw:
        graphs, gxy = df.calib_lam1n(zz_12ns[m-1], zz_1ns[m], lam1n, roi)
        if btn_detail.is_on():
            pf.graph_many(graphs, 'calibrate', gxy)

        zz_12ns[m], graphs, gxy = df.stitch(zz_12ns[m-1], zz_1ns[m], lam12, lam1n, roi)
        if btn_detail.is_on():
            pf.graph_many(graphs, 'stitch', gxy, sxy=(1, .75))

    if 2 <= m <= nw:
        zz_12ns[m] = df.cyclic_medfilter(zz_12ns[m], mnf, lam12)
        zz_proc = zz_12ns[m].copy()
        _, noise[m] = df.roi_measure(zz_12ns[m], roi)

    ent_np.set_entry(m)
    plot_n()

    print(f'> make_zz: m = {m}; sign = {sign}; wl_{m} = {wl:.8f}; lam_1{m} = {lam1n:.1f}; ' 
          f'zh = {zh:.1f}; z0_roi = {z0_roi:.1f}; mnf = {mnf}')


def plot_n():
    np = ent_np.get_val()
    nw = ent_nw.get_val()
    wl = wln[np]
    lam1n = lam_1ns[np]
    lam12 = lam_1ns[2]

    if np > nw: return

    if np == 0: rep = '_aa.png'; ulim = (0., 1.)
    else: rep = f'_{np}p.png'; ulim = pi2limit
    capA, _ = gf.path_parts(txt_path.replace('.txt', rep))
    pf.plotAAB(hhh[np], figname='HHn', capA=capA, roi=roi, sxy=sxy, ulimit=ulim)

    if np < 1: return

    pf.plotAAB(zz_ns[np], figname='ZZn', capA=f'ZZ_{np}', capB=f'lam_{np} = {wl:.8f}',
               roi=roi, sxy=sxy, ulimit=(-wl/2, wl/2))

    if np < 2: return

    pf.plotAAB(zz_1ns[np], figname='ZZ1n', capA=f'ZZ_1{np}', capB=f'lam_1{np} = {lam1n:.1f}',
               roi=roi, sxy=sxy, ulimit=(-lam1n/2, lam1n/2))

    pf.plotAAB(zz_12ns[np], figname='ZZ12n', capA=f'ZZ_12{np}',
               capB=f'lam12 = {lam12:.1f}; lam1{np} = {lam1n:.1f}; noise = {noise[np]:.1f}',
               roi=roi, sxy=sxy, ulimit=(-lam12/2, lam12/2))


def proc_zz():
    global zz_proc, proc_noise

    m = ent_n.get_val()
    lam12 = lam_1ns[2]
    lam1n = lam_1ns[m]
    # qxys = ent_qxy.get_list_val(float)
    # mnfp = ent_mnfp.get_list_val()

    zz_proc = df.zz_tilt(zz_12ns[m], nxydw, qxys, lam12)
    zz_proc = df.cyclic_medfilter(zz_proc, mnfp, lam12)

    # z0_roi = ent_z0roi.get_val(float)
    z_roi, _ = df.roi_cyclic_measure(zz_proc, roi, lam12)
    zz_proc = np.mod(zz_proc - z_roi + z0_roi + lam12/2, lam12) - lam12/2

    _, proc_noise = df.roi_measure(zz_proc, roi)
    capB = f'lam12 = {lam12:.1f}; lam1{m} = {lam1n:.1f}; noise = {proc_noise:.1f}'
    pf.plotAAB(zz_proc, figname='ZZproc', capA=f'ZZ_proc', capB=capB, roi=roi, sxy=sxy, ulimit=(-lam12/2, lam12/2))

    print(f'\n> proc_zz: qxys = {qxys}; mnfp = {mnfp}; proc_noise = {proc_noise:.1f}')


def denoise():
    global zz_proc

    gamma = ent_gamma.get_val(float)
    d_noise = proc_noise * gamma
    zz_proc = np.round(zz_proc / d_noise) * d_noise
    _, noise_new = df.roi_measure(zz_proc, roi)
    pf.plotAAB(zz_proc, figname='ZZproc', capA=f'ZZ_proc', capB=f'lam12 = {lam_1ns[2]:.1f}; noise = {noise_new:.1f}',
               roi=roi, sxy=sxy, ulimit=(-lam_1ns[2]/2, lam_1ns[2]/2))
    print(f'> denoise: d_noise = {d_noise:.1f}: new_noise = {noise_new:.1f}')


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
    # notes += '\n' + gf.prn_list('roi', roi, 0) + f'; z0_roi = {ent_z0roi.get_val(float):.1f}; zh = {zh:.1f}; ' \
    notes += '\n' + gf.prn_list('roi', roi, 0) + f'; z0_roi = {z0_roi:.1f}; zh = {zh:.1f}; ' \
                     + gf.prn_list('mnf', mnf, 0)[2:]
    # mnfp = ent_mnfp.get_list_val()
    notes += '\n' + gf.prn_list('qxys', qxys, 6) + gf.prn_list('mnfp', mnfp, 0)[2:]
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


def clear_mem():
    import gc
    global hhh, zz_ns, zz_1ns, zz_12ns

    del hhh, zz_ns, zz_1ns, zz_12ns
    gc.collect()
    print(f'> clreared: hhh, zz_ns, zz_1ns, zz_12ns')


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
btn_graphall.command(graph_all)
btn_mayavi.command(mayavi)
btn_plot.command(plot_n)
btn_clear.command(clear_mem)

btn_proc.command(proc_zz)
btn_dnoise.command(denoise)

btn_adios.command(adios)

btn_inv.on()

tloop = 10
fp.after(tloop, program_loop)
fp.mainloop()






