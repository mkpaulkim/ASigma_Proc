import numpy as np
# import pycubelib.plotting_functions as pf
# import pycubelib.general_functions as gf

pi = np.pi
pi2 = pi * 2


def renormalize(aa, alimit, blimit, btype=None):

    atype = type(aa)
    aa = aa.astype(float)
    if alimit:
        aa = np.clip(aa, alimit[0], alimit[1])
    else:
        alimit = (np.min(aa), np.max(aa))
    amin, amax = alimit
    bmin, bmax = blimit
    bb = bmin + (bmax - bmin) * (aa - amin) / (amax - amin)
    if btype:
        bb = bb.astype(btype)
    else:
        bb = bb.astype(atype)

    return bb


def calib_lam1n(zz12n, zz1n, lam1n, roi):

    nbin = 100

    ix, iy, rx, ry = roi
    dzz = np.mod(zz12n - zz1n + lam1n/2, lam1n) - lam1n/2
    dlam = lam1n / nbin
    histo, uu = np.histogram(dzz, bins=nbin, range=(-lam1n/2, lam1n/2))
    graphs = [(dzz[iy, :], (0, 0), 'dZZ', (0, len(dzz)), ()),
              (histo, (0, 1), 'histogram', (-lam1n/2, lam1n/2), ())]
    # pf.graph_many(graphs, 'calibrate lam1n', col_row=(1, 2), pause=1)
    gxy = (1, 2)

    return graphs, gxy


def stitch(zz12n_in, zz1n, lam12, lam1n, roi):

    yy = np.round(zz12n_in / lam1n) * lam1n
    zz = yy + zz1n
    dzz = zz - zz12n_in
    ezz = (np.abs(dzz) > (0.5 * lam1n)) * lam1n * np.sign(dzz)
    zz12n_out = np.mod(zz - ezz + lam12/2, lam12) - lam12/2

    ix, iy, rx, ry = roi
    ylimit = (-1.25 * lam12/2, 1.25 * lam12/2)
    graphs = []
    graphs += [(zz12n_in[iy, :], (0, 0), 'zz12n_in', (), ylimit)]
    graphs += [(zz1n[iy, :], (0, 1), 'zz1n', (), ylimit)]
    graphs += [(yy[iy, :], (0, 2), 'yy', (), ylimit)]
    graphs += [(zz[iy, :], (0, 3), 'zz', (), ylimit)]
    graphs += [(dzz[iy, :], (0, 4), 'dzz', (), ylimit)]
    graphs += [(ezz[iy, :], (0, 5), 'ezz', (), ylimit)]
    graphs += [(zz12n_out[iy, :], (0, 6), 'zz12n_out', (0, len(zz)), ylimit)]
    # pf.graph_many(graphs, 'stitch', (1, 7), sxy=(.25, .25), pause=1)
    gxy = (1, 7)

    return zz12n_out, graphs, gxy


def diffract(hhp, hha, wl, nxydw, qxysz):
    nx, ny, dx, nw = nxydw
    qx, qy, qs, zh = qxysz
    kk = pi2 / wl
    ax, ay = (nx * dx, ny * dx)
    x = np.linspace(-ax/2, ax/2 - dx, nx)
    y = np.linspace(-ay/2, ay/2 - dx, ny)
    xx, yy = np.meshgrid(x, y)
    ak = pi2 / dx
    dkx, dky = (pi2 / ax, pi2 / ay)
    kx = np.linspace(-ak/2, ak/2 - dkx, nx)
    ky = np.linspace(-ak/2, ak/2 - dky, ny)
    kxx, kyy = np.meshgrid(kx, ky)

    hh = hha * np.exp(1j * hhp)
    ggq = np.exp(1j * kk * (xx * np.sin(qx) + yy * np.sin(qy) + (xx**2 + yy**2) * (qs / 2)))
    ggk = np.exp(1j * zh * np.sqrt(kk**2 - kxx**2 - kyy*2))
    ff = np.fft.fftshift(np.fft.fft2(hh * ggq))
    hh_out = np.fft.ifft2(np.fft.ifftshift(ff * ggk))
    hhp_out = np.angle(hh_out)
    hha_out = np.abs(hh_out)

    return hhp_out, hha_out


# def calib_lam1n_0(zz12n_in, zz1n, lam12, lam1n, roi):
#     nbin = 100
#     dlam = lam12 / nbin
#     ix, iy, rx, ry = roi
#     ep1n = zz1n * pi2 / lam1n
#     lam12limit = (-lam12/2, lam12/2)
#
#     print(f'>>> lam1n_in = {lam1n:.1f}')
#     while True:
#         ans = float(input('> lam1n = '))
#         if ans:
#             lam1n = ans
#         else:
#             break
#
#         zp1n = ep1n * lam1n / pi2
#         dzz = zz12n_in - zp1n
#         dzz1n = np.mod(dzz + lam1n/2, lam1n) - lam1n/2
#         # pf.plotAAB(dzz, capA=f'dzz: lam1n = {lam1n:.1f}', sxy=(.35, .35), pause=1)
#
#         # histo, uu = np.histogram(dzz, bins=nbin, range=(-lam12/2, lam12/2))
#         histo, uu = np.histogram(dzz1n, bins=nbin, range=(-lam1n/2, lam1n/2))
#
#         graphs = [(zz12n_in[iy, :], (0, 0), f'zz12n: lam12 = {lam12:.1f}', lam12limit),
#                   (zp1n[iy, :], (0, 1), f'zz1n_: lam1n = {lam1n:.1f}', lam12limit),
#                   (dzz[iy, :], (0, 2), f'dzz', lam12limit),
#                   (dzz1n[iy, :], (0, 3), f'dzz1n', (-lam1n/2, lam1n/2))]
#         pf.graph_many(graphs, 'calib', col_row=(1, 4), sxy=(.35, .3), pause=1)
#         # pf.graphB(histo, caption='histogram', xpars=(-lam12/2, lam12/nbin), sxy=(7.5, .3), line='-+', pause=1)
#         pf.graphB(histo, caption='histogram', xpars=(-lam1n/2, lam1n/nbin), sxy=(7., .3), line='-+', pause=1)
#
#         zz_12n = stitch(zz12n_in, zz1n, lam12, lam1n)
#         pf.plotAAB(zz_12n, capA=f'ZZ12n: lam_1n = {lam1n:.1f}', roi=roi, sxy=(.35, .35), pause=1)
#
#         # xx = np.arange(nbin) * lam12 / nbin - lam12/2
#         # pf.plt.figure('histogram')
#         # pf.plt.plot(xx, histo, '-+')
#         # pf.plt.pause(1)
#
#     return lam1n


def cyclic_medfilter(zz, mnf, lamz):
    import scipy.signal as sig

    mf, nf = mnf
    uu = zz * pi2 / lamz
    cc = np.cos(uu)
    ss = np.sin(uu)
    for n in range(nf):
        cc = sig.medfilt2d(cc, mf)
        ss = sig.medfilt2d(ss, mf)
    zz_out = np.arctan2(ss, cc) * lamz / pi2

    return zz_out


if __name__ == '__main__':

    pass




