import numpy as np
import pycubelib.tkinter_parts as tp
import pycubelib.files_functions as ff
import pycubelib.general_functions as gf
import pycubelib.plotting_functions as pf
import pycubelib.data_functions as df

pi = np.pi
pi2 = pi * 2
pi2limit = (-pi, pi)
sxy = (.85, 1.)
view = (20, -10)

fp = tp.tkwindow('GammaSigmaProc', (2030*1 + 50, 50, 1000, 400), tkbg='gray85')
btn_readtxt = tp.CmdButton(fp, (50, 50, 10), 'read TXT')
ent_txtpath = tp.ParamEntry(fp, (170, 55, 35), '', '')
ent_nxydw = tp.ParamEntry(fp, (170, 100, 35), 'nxydw', '')
ent_lams = tp.ParamEntry(fp, (170, 150, 35), 'lam12, lam1n', '0.0, 0.0')
btn_readhhp = tp.CmdButton(fp, (50, 200, 10), 'read HHp')
ent_hhppath = tp.ParamEntry(fp, (170, 200, 35), '', '')
ent_n = tp.ParamEntry(fp, (100, 250, 5), 'n', '0')
prog_n = tp.ProgressBar(fp, (170, 250, 285), '')
ent_roi = tp.ParamEntry(fp, (170, 300, 25), 'ixy, rxy', '200, 550, 30, 30')
ent_z0roi = tp.ParamEntry(fp, (550, 50, 15), 'roi_Z0', '-1500.0')
ent_mnf = tp.ParamEntry(fp, (550, 100, 15), 'mf, nf', '3, 1')
ent_qxy = tp.ParamEntry(fp, (550, 150, 25), 'qxys', '0.000, 0.000, 0.000') # qx, qy in rad; ss curvature in 1/mm
ent_pmnf = tp.ParamEntry(fp, (550, 200, 15), 'p_mf, nf', '3, 1')
btn_betascan = tp.CmdButton(fp, (800, 50, 10), 'beta scan')
btn_gammascan = tp.CmdButton(fp, (800, 100, 10), 'gamma scan')
btn_proc = tp.CmdButton(fp, (800, 150, 10), 'proc')
btn_mayavi = tp.CmdButton(fp, (800, 200, 10), 'mayavi')
btn_inv = tp.CmdButton(fp, (550, 250, 5), 'inv', 'gray90')
# ent_ph = tp.ParamEntry(fp, (550, 300, 15), 'ph_pi', '0.0')

btn_adios = tp.CmdButton(fp, (800, 300, 10), 'adios', 'indian red')

txt_path = ''
nxydw = []
wln = []
lam_1ns = [0, 0, 0]
lam12 = 0
blank = []
hhh = []
roi = []
z0_roi = 0
mnf = []
qxys = []
pmnf = []
zz = []
zz_proc = []
phis_roi = []


def program_loop():
    global roi, z0_roi, mnf, qxys, pmnf, lam12

    roi = ent_roi.get_list_val()
    z0_roi = ent_z0roi.get_val(float)
    mnf = ent_mnf.get_list_val()
    qxys = ent_qxy.get_list_val(float)
    pmnf = ent_pmnf.get_list_val()
    lam12 = lam_1ns[2]

    fp.after(tloop, program_loop)


def read_txt():
    global txt_path, nxydw, wln

    text, txt_path = ff.read_txt()
    t = text.find('%%%')
    if t < 0:
        t = text.find('>>>')
    param_txt = text[t:]
    param_txt.replace('([', '[')            # for 2018-01 Saturn data
    nx = gf.find_param(param_txt, 'nx')
    ny = gf.find_param(param_txt, 'ny')
    dx = gf.find_param(param_txt, 'dx', float)
    nw = gf.find_param(param_txt, 'nw')
    wln = gf.find_param(param_txt, 'wln', float)
    nxydw = (nx, ny, dx, nw)

    ent_txtpath.set_entry(txt_path)
    ent_nxydw.set_entry(nxydw)

    txt_note = f'>> txt_path = {txt_path}'
    txt_note += f'\n{text}'
    if len(wln) < nw + 1:
        wln = [0.0] + wln
        txt_note += f'\n> !!! len(wln) < nw + 1 = {nw + 1}'
        txt_note += f'\n> new wln: {gf.prn_list("wl", wln)} \n'

    print(txt_note)
    print(f'> read TXT: done ...')

    initialize()


def initialize():
    global blank, hhh, lam_1ns, zz

    ent_hhppath.set_entry('')
    ent_n.set_entry(0)
    prog_n.setval(0)

    nx, ny, dx, nw = nxydw
    blank = np.zeros((ny, nx))
    zz = np.copy(blank)
    hhh = []

    lam_1ns = [0, 0]
    for n in range(2, nw+1):
        lam_1ns += [wln[1] * wln[n] / (wln[1] - wln[n])]

    ent_lams.set_entry(f'{lam_1ns[2]:.1f}, {lam_1ns[-1]:.1f}')


def read_hhp():
    global hhh
    nx, ny, dx, nw = nxydw

    hhh = []
    for n in range(0, nw+1):
        if n == 0:
            hh_path = txt_path.replace('.txt', f'_aa.png')
            hh, _ = ff.read_png(hh_path, alimit=(0., 1.))
            fig_name = 'figure 1'
        else:
            hh_path = txt_path.replace('.txt', f'_{n}p.png')
            hh, _ = ff.read_png(hh_path, alimit=pi2limit)
            fig_name = 'figure 2'
        hhh += [hh]

        ent_hhppath.set_entry(hh_path)
        ent_n.set_entry(n)
        prog_n.setval(100 * n / nw)

        capA, _ = gf.path_parts(hh_path)
        pf.plotAAB(hh, figname=fig_name, capA=capA, roi=roi, sxy=sxy)

        print(f'> hh_path = {hh_path}')

    print(f'> read HHp: done ...')


def beta_scan():
    global zz, phis_roi

    nx, ny, dx, nw = nxydw
    sumphi = np.copy(blank)
    limit_12 = (-lam12/2, lam12/2)
    phis_roi = [0, 0]

    for n in range(2, nw+1):
        ent_n.set_entry(n)
        prog_n.setval(100 * n / nw)
        lam1n = lam_1ns[n]
        lamnn = wln[n - 1] * wln[n] / (wln[n - 1] - wln[n])

        delphi = np.mod(hhh[n-1] - hhh[n] + pi, pi2) - pi
        delphi = df.cyclic_medfilter(delphi, mnf, pi2)
        ph0_roi = z0_roi * pi2 / lam12

        ph_roi, _ = df.roi_cyclic_measure(delphi, roi, pi2)
        delphi = np.mod(delphi - ph_roi + ph0_roi + pi, pi2) - pi
        _, del_noise = df.roi_measure(delphi * lamnn / pi2, roi)
        # pf.plotAAB(delphi * lamnn / pi2, figname='figure 1', capA=f'Z_delphi: n = {n}', roi=roi,
        #            sxy=sxy, ulimit=limit_12, capB=f'lamnn = {lamnn:.1f}; noise = {del_noise:.1f}')
        phis_roi += [ph_roi - ph0_roi]
        # phis_roi += [ph_roi]

        sumphi += delphi
        zz = sumphi * lam1n / pi2
        _, sum_noise = df.roi_measure(zz, roi)
        # pf.plotAAB(zz, figname='figure 2', capA=f'ZZ: beta_sacn: n = {n}', roi=roi,
        #            sxy=sxy, ulimit=limit_12, capB=f'lam1n = {lam1n:.1f}; noise = {sum_noise:.1f}')

        print(f'> ZZ: beta_scan: n = {n}: del_noise = {del_noise:.1f}; sum_noise = {sum_noise:.1f}')


def gamma_scan():
    global zz

    nx, ny, dx, nw = nxydw
    limit_12 = (-lam12/2, lam12/2)
    lam1n = lam_1ns[-1]
    eee = np.zeros((ny, nx, nw)) * 1j
    hha = hhh[0]
    phi = 0
    # phi = ent_ph.get_val(float) * pi

    for n in range(1, nw+1):
        ent_n.set_entry(n)
        prog_n.setval(100 * n / nw)
        phi += phis_roi[n]

        eee[:, :, n-1] = hha * np.exp(1j * hhh[n]) * np.exp(1j * phi)

    fff = np.fft.fft(eee, axis=2)
    fff = np.fft.fftshift(fff, axes=2)
    fffa = np.abs(fff)
    pff = np.transpose(np.log10(fffa[roi[1], :, :]))

    z_step = lam12 / (nw-1)
    zz_ffa = - 1.0 * (np.argmax(fffa, axis=2) * z_step - lam12/2)
    # zz_ffa = df.cyclic_medfilter(zz_ffa, mnf, lam12)
    # z_roi, _ = df.roi_cyclic_measure(zz_ffa, roi, lam12)
    # zz_ffa = np.mod(zz_ffa - z_roi + z0_roi + lam12/2, lam12) - lam12/2

    zz1n = (hhh[1] - hhh[-1]) * lam1n / pi2
    zz1n = df.cyclic_medfilter(zz1n, mnf, lam1n)
    z_1n, _ = df.roi_cyclic_measure(zz1n, roi, lam1n)
    zz1n = np.mod(zz1n - z_1n + z0_roi + lam1n/2, lam1n) - lam1n/2
    _, noise1n = df.roi_measure(zz1n, roi)

    zz = np.mod(zz_ffa + zz1n + lam12/2, lam12) - lam12/2
    zz = df.cyclic_medfilter(zz, mnf, lam12)
    _, noise = df.roi_measure(zz, roi)

    pf.plotAAB(pff, figname='figure 1', capA='log_ffa', sxy=sxy, cmap='jet', crsr=False)
    pf.plotAAB(zz_ffa, figname='figure 2', capA='ZZ_ffa', capB=f' ', ulimit=limit_12, sxy=sxy, roi=roi)
    pf.plotAAB(zz1n, figname='figure 3', capA='ZZ_1n', capB=f'noise_1n = {noise1n:.1f}', ulimit=(-lam1n/2, lam1n/2),
               sxy=sxy, roi=roi)
    pf.plotAAB(zz, figname='figure 4', capA='ZZ: gamma_scan', capB=f'noise = {noise:.1f}', ulimit=limit_12,
               sxy=sxy, roi=roi)

    print(f'> ZZ: gamma_scan: noise = {noise:.1f}')


def proc_zz():
    global zz_proc
    limit_12 = (-lam12/2, lam12/2)

    zz_proc = df.zz_tilt(zz, nxydw, qxys, lam12)
    zz_proc = df.cyclic_medfilter(zz_proc, pmnf, lam12)
    z_roi, _ = df.roi_cyclic_measure(zz_proc, roi, lam12)
    zz_proc = np.mod(zz_proc - z_roi + z0_roi + lam12/2, lam12) - lam12/2
    if btn_inv.is_on():
        zz_proc = - 1. * zz_proc

    _, proc_noise = df.roi_measure(zz_proc, roi)
    capB = f'noise = {proc_noise:.1f}'
    pf.plotAAB(zz_proc, figname='figure 5', capA=f'ZZ_proc', capB=capB, roi=roi, sxy=sxy, ulimit=limit_12, cmap='jet')

    print(f'\n> ZZ_proc: qxys = {qxys}; pmnf = {pmnf}; proc_noise = {proc_noise:.1f}')


def mayavi():
    cap = gf.path_parts(txt_path)[0]
    pf.mayaviAA(zz_proc, caption=cap, ulimit=(-lam12/2, lam12/2), view=view, cmap='jet')

    print(f'> mayavi: done ...')


def adios():
    fp.destroy()
    print(f'> adios amigos ...')
    quit()


btn_readtxt.command(read_txt)
btn_readhhp.command(read_hhp)
btn_betascan.command(beta_scan)
btn_gammascan.command(gamma_scan)
btn_proc.command(proc_zz)
btn_inv.command(btn_inv.switch)
btn_mayavi.command(mayavi)

btn_adios.command(adios)

tloop = 10
fp.after(tloop, program_loop)
fp.mainloop()


'''

sxy = (.85, 1)
view = (20, -10)

fp = tp.tkwindow('GammaSigmaProc', (2030*0 + 50, 50, 1400, 400), tkbg='gray85')
btn_readtxt = tp.CmdButton(fp, (50, 50, 10), 'read TXT')
ent_txtpath = tp.ParamEntry(fp, (170, 55, 35), '', '')
btn_readhhp = tp.CmdButton(fp, (50, 100, 10), 'read HHp')
ent_hhppath = tp.ParamEntry(fp, (170, 105, 35), '', '')

btn_qfft = tp.CmdButton(fp, (50, 200, 10), 'do qfft')
btn_note = tp.CmdButton(fp, (850, 150, 10), 'notes')
btn_save = tp.CmdButton(fp, (850, 200, 10), 'save')
btn_adios = tp.CmdButton(fp, (850, 50, 10), 'adios', 'indian red')

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
fzz1n = []


def program_loop():
    global roi, z0_roi, zh, mnf, qxys, mnfp

    roi = ent_roi.get_list_val()
    z0_roi = ent_z0roi.get_val(float)
    zh = ent_zh.get_val(float)
    mnf = ent_mnf.get_list_val()
    qxys = ent_qxy.get_list_val(float)
    mnfp = ent_mnfp.get_list_val()

    fp.after(tloop, program_loop)


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
zz
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


btn_qfft.command(do_qfft)
btn_note.command(print_notes)
btn_save.command(save_zzproc)
btn_inv.command(btn_inv.switch)
btn_mayavi.command(mayavi)
btn_proc.command(proc_zz)



'''



