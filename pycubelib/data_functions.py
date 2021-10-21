import numpy as np
import pycubelib.plotting_functions as pf


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


def get_roi(aa, roi):
    rx0, ry0, rx, ry = roi
    return aa


def stitch(zz12n_in, zz1n, lam12, lam1n):
    yy = np.round(zz12n_in / lam1n) * lam1n
    zz = yy + zz1n
    dzz = zz - zz12n_in
    ezz = (np.abs(dzz) > (0.5 * lam1n)) * lam1n * np.sign(dzz)
    zz12n_out = np.mod(zz - ezz + lam12/2, lam12) - lam12/2

    if 1:
        iy = 950
        ylimit = (-1.5 * lam12/2, 1.5 * lam12/2)
        graphs = []
        graphs += [(zz12n_in[iy, :], (0, 0), 'zz12n_in', ylimit)]
        graphs += [(zz1n[iy, :], (0, 1), 'zz1n', ylimit)]
        graphs += [(yy[iy, :], (0, 2), 'yy', ylimit)]
        graphs += [(zz[iy, :], (0, 3), 'zz', ylimit)]
        graphs += [(dzz[iy, :], (0, 4), 'dzz', ylimit)]
        graphs += [(ezz[iy, :], (0, 5), 'ezz', ylimit)]
        graphs += [(zz12n_out[iy, :], (0, 6), 'zz12n_out', ylimit)]
        pf.graph_many(graphs, 'stitch', (1, 7), sxy=(.25, .25), pause=1)

    return zz12n_out


def stitch_x(zz12n_in, zz1n, lam12, lam1n):
    yy = np.round(zz12n_in / lam1n) * lam1n
    zz = np.mod(yy + zz1n + lam12/2, lam12) - lam12/2
    dzz = np.mod(zz - zz12n_in + lam12/2, lam12) - lam12/2
    ezz = (np.abs(dzz) > (0.5 * lam1n)) * lam1n * np.sign(dzz)
    zz12n_out = np.mod(zz - ezz + lam12/2, lam12) - lam12/2
    return zz12n_out


if __name__ == '__main__':

    pass




