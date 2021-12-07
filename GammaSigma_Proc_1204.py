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
view = (20, -10)

fp = tp.tkwindow('BetaSigmaProc', (20, 50, 1400, 400), tkbg='gray85')

btn_readtxt = tp.CmdButton(fp, (100, 50, 10), 'read TXT')
btn_readphs = tp.CmdButton(fp, (100, 100, 10), 'read HHp')
btn_qfft = tp.CmdButton(fp, (100, 200, 10), 'do qfft')
btn_note = tp.CmdButton(fp, (850, 150, 10), 'notes')
btn_save = tp.CmdButton(fp, (850, 200, 10), 'save')
btn_adios = tp.CmdButton(fp, (850, 50, 10), 'adios', 'indian red')

ent_txtpath = tp.ParamEntry(fp, (220, 55, 35), '', '')
ent_phspath = tp.ParamEntry(fp, (220, 105, 35), '', '')
ent_n = tp.ParamEntry(fp, (250, 205, 5), 0, 'n')
ent_nw = tp.ParamEntry(fp, (350, 205, 5), 0, 'nw')
ent_lam1n = tp.ParamEntry(fp, (300, 250, 15), 0, 'lam1n')
ent_roi = tp.ParamEntry(fp, (600, 50, 25), '200, 750, 10, 10', 'ixy rxy')
ent_z0roi = tp.ParamEntry(fp, (600, 100, 15), 0.0, 'Z0_roi')
ent_zh = tp.ParamEntry(fp, (600, 150, 15), 0.0, 'Zh')
ent_mnf = tp.ParamEntry(fp, (600, 250, 15), '3, 1', 'mnf')

prog_n = tp.ProgressBar(fp, (100, 310, 410), '')

btn_proc = tp.CmdButton(fp, (1000, 50, 10), 'proc')
ent_qxy = tp.ParamEntry(fp, (1000, 100, 25), '0.000, 0.000, 0.000', 'qxys') # qx, qy in rad; ss curvature in 1/mm
ent_mnfp = tp.ParamEntry(fp, (1000, 150, 15), '3, 10', 'mnfp')
btn_inv = tp.CmdButton(fp, (1300, 150, 5), 'inv', 'gray90')
btn_mayavi = tp.CmdButton(fp, (1000, 300, 10), 'mayavi')

txt_path = ''
roi = (0, 0, 0, 0)
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
fzz = []
fzz1n=[]


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
    if t < 0:
        t = text.find('>>>')
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
    print(f'> read TXT: done ...')

    get_lam1ns()


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

    print('\n>>> get_lam1ns: reset zz_ns, zz_1ns, zz_12ns ...')
    print(gf.prn_list('lam_1ns', lam_1ns, 1) + '\n')


def do_qfft():
    global zz_ns, zz_1ns, zz_12ns, zz_proc, noise, fzz, fzz1n

    nx, ny, dx, nw = nxydw
    # phi = np.copy(blank)
    lam12 = lam_1ns[2]
    lam1n = lam_1ns[nw]
    eee = np.zeros((ny, nx, nw+1)) * 1j
    hha = hhh[0]
    phi = 0.

    for n in range(1, nw+1):
        ent_n.set_entry(n)
        # ent_lam1n.set_entry(lam_1ns[n+1])
        print(f'> set_phase: hhh[{n}]')
        ent_n.set_entry(n)
        prog_n.setval(100 * n / nw)
        ph0_roi = np.mod(z0_roi * pi2 / lam12 + pi, pi2) - pi

        if n == 1:
            eee[:, :, 1] = hha * np.exp(1j * hhh[1])
        else:
            delphi = np.mod(hhh[n] - hhh[n-1] + pi, pi2) - pi
            # delphi = df.cyclic_medfilter(delphi, mnf, pi2)
            ph_roi, _ = df.roi_cyclic_measure(delphi, roi, pi2)
            delphi = np.mod(delphi - ph_roi + ph0_roi + pi, pi2) - pi
            pf.plotAAB(delphi, figname='delphi', capA=f'hhh[{n}] - hhh[{n+1}]', roi=roi, sxy=sxy, ulimit=(-pi, pi))

            phi += ph_roi - ph0_roi
            # phi = 0
            eee[:, :, n] = hha * np.exp(1j * hhh[n]) * np.exp(- 1j * phi)

    eee = eee[:, :, 1:]
    fff = np.fft.fft(eee, axis=2)
    fff = np.fft.fftshift(fff, axes=2)
    fffa = np.abs(fff)
    pf.plotAAB(fffa[:, roi[0], :], figname='fffa', roi=roi, sxy=(.85, nw/2000), crsr=False)

    dz = lam12 / (nw-1)
    fzz = np.argmax(fffa, axis=2) * dz - lam12/2
    pf.plotAAB(fzz, figname='fzz', roi=roi, sxy=sxy)

    # zz1n = (np.mod(hhh[nw] - hhh[1] + pi, pi2) - pi) * lam1n / pi2
    zz1n = (hhh[nw] - hhh[1]) * lam1n / pi2
    # zz1n = df.cyclic_medfilter(zz1n, mnf, lam1n)
    pf.plotAAB(zz1n, figname='zz1n', roi=roi, sxy=sxy)

    fzz_roi, _ = df.roi_cyclic_measure(fzz, roi, lam12)
    fzz1n_roi = np.mod(fzz_roi + lam1n/2, lam1n) - lam1n/2
    zz1n_roi, _ = df.roi_cyclic_measure(zz1n, roi, lam1n)
    fzz1n = fzz + zz1n + fzz1n_roi - zz1n_roi
    pf.plotAAB(fzz1n, figname='fzz1n', roi=roi, sxy=sxy)



def proc_zz():
    global zz_proc, proc_noise

    m = ent_n.get_val()
    lam12 = lam_1ns[2]
    lam1n = lam_1ns[m]
    zz_proc = np.copy(fzz1n)

    zz_proc = df.zz_tilt(zz_proc, nxydw, qxys, lam12)
    zz_proc = df.cyclic_medfilter(zz_proc, mnfp, lam12)

    z_roi, _ = df.roi_cyclic_measure(zz_proc, roi, lam12)
    zz_proc = np.mod(zz_proc - z_roi + z0_roi + lam12/2, lam12) - lam12/2

    if btn_inv.is_on():
        zz_proc = - 1. * zz_proc

    _, proc_noise = df.roi_measure(zz_proc, roi)
    capB = f'lam12 = {lam12:.1f}; lam1{m} = {lam1n:.1f}; noise = {proc_noise:.1f}'
    pf.plotAAB(zz_proc, figname='ZZproc', capA=f'ZZ12_{m}_proc', capB=capB, roi=roi, sxy=sxy, ulimit=(-lam12/2, lam12/2))

    print(f'\n> proc_zz: qxys = {qxys}; mnfp = {mnfp}; proc_noise = {proc_noise:.1f}')


def mayavi():
    n = ent_n.get_val()
    cap = gf.path_parts(txt_path)[0] + f': ZZ12_{n}_proc'
    lam12 = lam_1ns[2]
    pf.mayaviAA(zz_proc, caption=cap, ulimit=(-lam12/2, lam12/2), view=view)

    print(f'> mayavi: done ...')


def print_notes():
    global notes
    pass

    # notes = txt_note
    # notes += f'> ================'
    # notes += '\n' + gf.runstamp(__file__)
    # notes += '\n' + gf.prn_list('lam_1ns', lam_1ns, 1)
    # # notes += '\n' + gf.prn_list('roi', roi, 0) + f'; z0_roi = {ent_z0roi.get_val(float):.1f}; zh = {zh:.1f}; ' \
    # notes += '\n' + gf.prn_list('roi', roi, 0) + f'; z0_roi = {z0_roi:.1f}; zh = {zh:.1f}; ' \
    #                  + gf.prn_list('mnf', mnf, 0)[2:]
    # # mnfp = ent_mnfp.get_list_val()
    # notes += '\n' + gf.prn_list('qxys', qxys, 6) + gf.prn_list('mnfp', mnfp, 0)[2:]
    # notes += '\n' + gf.prn_list('noise', noise, 1)[:-2] + f' + [{proc_noise:.1f}]; '
    #
    # print()
    # print(f'> notes ---------------------------------------------------------------------')
    # print(notes)
    # print(f'> ---------------------------------------------------------------------------')
    # print()


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
btn_qfft.command(do_qfft)
btn_note.command(print_notes)
btn_save.command(save_zzproc)
btn_inv.command(btn_inv.switch)
btn_mayavi.command(mayavi)
btn_proc.command(proc_zz)

btn_adios.command(adios)

tloop = 10
fp.after(tloop, program_loop)
fp.mainloop()






