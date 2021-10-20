import numpy as np
import pycubelib.general_functions as gf


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
    zz = np.mod(yy + zz1n + lam12/2, lam12) - lam12/2
    dzz = np.mod(zz - zz12n_in + lam12/2, lam12) - lam12/2
    ezz = (np.abs(dzz) > (0.5 * lam1n)) * lam1n * np.sign(dzz)
    zz12n_out = np.mod(zz - ezz + lam12/2, lam12) - lam12/2
    return zz12n_out


if __name__ == '__main__':

    pass




