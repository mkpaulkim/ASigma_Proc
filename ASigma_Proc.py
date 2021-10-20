import numpy as np
from tkinter import filedialog
import pycubelib.files_functions as ff
import pycubelib.general_functions as gf
import pycubelib.plotting_functions as pf

pi = np.pi
pi2 = pi * 2
pi2limit = (-pi, pi)

txt_path = filedialog.askopenfilename(title='TXT file path', filetypes=[('txt files', '*.txt')])
print(f'> txt_path = {txt_path}')
notes, txt_path = ff.read_txt(txt_path)
print(f'> notes: \n{notes}')
nx = gf.find_param(notes, 'nx', int)
ny = gf.find_param(notes, 'ny', int)
dx = gf.find_param(notes, 'dx_um', float)
nh = gf.find_param(notes, 'nh', int)
lam_ns = gf.find_param(notes, 'lam_ns', float)

sxy = (.3, .3)
hhhp = [np.zeros((ny, nx))]
for n in range(nh+1):
    if n == 0:
        hh_path = txt_path.replace('.txt', f'_aa.png')
    else:
        hh_path = txt_path.replace('.txt', f'_{n}p.png')
    hhp, _ = ff.read_png(hh_path, alimit=pi2limit)
    capA, _ = gf.path_parts(hh_path)
    pf.plotAAB(hhp, f'hh{n}p', capA=capA, sxy=sxy, pause=.5)
    hhhp += [hhp]
hha = hhhp[0]

zz12 = np.mod(hhhp[1] - hhhp[2] + pi, pi2) - pi
pf.plotAAB(zz12, 'zz12', capA='ZZ12', sxy=sxy)

lam_1ns = np.zeros(nh + 1)
lam1 = lam_ns[1]
for n in range(2, nh+1):
    lam_1ns[n] = (lam1 * lam_ns[n]) / (lam1 - lam_ns[n])
print(gf.prn_list('lam_1ns', lam_1ns, 1))

pf.plt.show()



'''
def program_loop():
    global g_t

    if fp.btn_setcam.is_on():
        m.set_cam()
        fp.btn_setcam.off()

    if fp.btn_camera.is_on():
        g_t, g_cc = m.show_cam(g_t)

    if fp.btn_angles.is_on():
        m.calc_angles()
        fp.btn_angles.off()

    if fp.btn_vset.is_on():
        v_piezo = fp.ent_vpiezo.get_val(float)
        hw.piezo.set_vpz(v_piezo)
        fp.btn_vset.off()

    if fp.btn_hpath.is_on():
        txt_path = filedialog.asksaveasfilename(title='HH.txt file path', filetypes=[('txt files', '*.txt')])
        # txt_path = filedialog.askopenfilename(title='HH.txt file path', filetypes=[('txt files', '*.txt')])
        if txt_path.find('.txt') < 0: txt_path += '.txt'
        fp.ent_hfile.set_entry(txt_path)
        if fp.btn_readh.is_on(): m.read_hparams()
        fp.btn_hpath.off()

    if fp.btn_bpath.is_on():
        txt_path = filedialog.askopenfilename(title='BB.txt file path', filetypes=[('txt files', '*.txt')])
        if txt_path.find('.txt') < 0: txt_path += '.txt'
        fp.ent_bfile.set_entry(txt_path)
        fp.btn_bpath.off()

    if fp.btn_scanh.is_on():
        m.start_scan()
        # fp.btn_scanh.off()

    if fp.btn_procz.is_on():
        m.proc_zzz()
        fp.btn_procz.off()

    if fp.btn_cfg.is_on():
        m.close_figs()
        fp.btn_cfg.off()

    if fp.btn_graphall.is_on():
        c = f'{fp.ent_cols.get_val()}'
        cc = []
        for i in range(len(c)):
            cc += [int(c[i])]
        cc = tuple(cc)
        nrow = fp.ent_nrow.get_val()
        pf.graph_many('graph_all', m.g_graphs, cols=cc, nrow=nrow, resize=0.5)
        fp.btn_graphall.off()

    fp.tkwindow.after(tloop, program_loop)


def adios():
    hw.qcam.dispose()
    fp.tkwindow.destroy()
    print(f'> adios amigos ...')
    quit()


fp.btn_adios.command(adios)
g_t = hw.qcam.tframe

tloop = 10
fp.tkwindow.after(tloop, program_loop)
fp.tkwindow.mainloop()

'''



