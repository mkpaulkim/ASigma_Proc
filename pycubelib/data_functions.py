import numpy as np

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
    histo, uu = np.histogram(dzz, bins=nbin, range=(-lam1n/2, lam1n/2))
    graphs = [(dzz[iy, :], (0, 0), 'dZZ', (0, len(dzz)), ()),
              (histo, (0, 1), 'histogram', (-lam1n/2, lam1n/2), ())]
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
    gxy = (1, 7)

    return zz12n_out, graphs, gxy


def diffract(hhp, hha, wl, nxydw, zh):
    nx, ny, dx, nw = nxydw
    kk = pi2 / wl
    ax, ay = (nx * dx, ny * dx)
    ak = pi2 / dx
    dkx, dky = (pi2 / ax, pi2 / ay)
    kx = np.linspace(-ak/2, ak/2 - dkx, nx)
    ky = np.linspace(-ak/2, ak/2 - dky, ny)
    kxx, kyy = np.meshgrid(kx, ky)

    hh = hha * np.exp(1j * hhp)
    ggk = np.exp(1j * zh * np.sqrt(kk**2 - kxx**2 - kyy*2))
    ff = np.fft.fftshift(np.fft.fft2(hh))
    hh_out = np.fft.ifft2(np.fft.ifftshift(ff * ggk))
    hhp_out = np.angle(hh_out)
    hha_out = np.abs(hh_out)

    return hhp_out, hha_out


def zz_tilt(zz, nxydw, qxys, zlam):
    nx, ny, dx, nw = nxydw
    qx, qy, qs = qxys
    qs *= 1e-3  # convert curvature from 1/mm to 1/um
    ax, ay = (nx * dx, ny * dx)
    x = np.linspace(-ax/2, ax/2 - dx, nx)
    y = np.linspace(-ay/2, ay/2 - dx, ny)
    xx, yy = np.meshgrid(x, y)
    zzq = (xx * np.sin(qx) + yy * np.sin(qy) + (xx**2 + yy**2) * (qs / 2))
    zz_out = np.mod(zz + zzq + zlam/2, zlam) - zlam/2
    return zz_out


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


def get_aaroi(aa, roi):
    ix, iy, rx, ry = roi
    rx1 = ix - rx//2
    rx2 = ix + rx//2
    ry1 = iy - ry//2
    ry2 = iy + ry//2
    aa_roi = aa[ry1:ry2, rx1:rx2]
    return aa_roi


def roi_measure(aa, roi):
    aa_roi = get_aaroi(aa, roi)
    roi_ave = np.mean(aa_roi)
    roi_std = np.std(aa_roi)
    return roi_ave, roi_std


def roi_cyclic_measure(zz, roi, zlambda):
    zz_roi = get_aaroi(zz, roi)
    cc = np.cos(zz_roi * pi2 / zlambda)
    ss = np.sin(zz_roi * pi2 / zlambda)
    c_ave = np.mean(cc)
    s_ave = np.mean(ss)
    cyc_ave = np.arctan2(s_ave, c_ave) * zlambda / pi2
    c_med = np.median(cc)
    s_med = np.median(ss)
    cyc_med = np.arctan2(s_med, c_med) * zlambda / pi2
    return cyc_ave, cyc_med


if __name__ == '__main__':

    pass




