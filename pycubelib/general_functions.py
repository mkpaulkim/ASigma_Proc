import numpy as np
import datetime
import os

pi = np.pi
pi2 = pi * 2


def runstamp(script_path):
    s_path, _ = path_parts(script_path, 2)
    time_stamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    stamp = f'[.../{s_path}] {time_stamp}'
    return stamp


def path_parts(full_path, nlast=2):
    base_path, fname_ext = os.path.split(full_path)
    fname, ext = os.path.splitext(fname_ext)

    rest_path = full_path
    parts = []
    for n in range(nlast):
        rest_path, part = os.path.split(rest_path)
        parts = [part] + parts
    short_path = f'{"/".join(parts)}'

    parts_dict = {'full_path': full_path, 'base_path': base_path, 'fname': fname, 'ext': ext, 'short_path': short_path}

    return short_path, parts_dict


def what_is(name, var, cmax=200):
    if type(var) == str:
        vlen = f' of len {len(var)}'
    elif type(var) == list:
        vlen = f' of len {len(var)}'
    elif type(var) == np.ndarray:
        vlen = f' of size {var.shape} of element type {type(var.reshape(np.product(var.shape))[0])}'
    else:
        vlen = ''

    output = f'> {name} is {type(var)}' + vlen + f': {name} = {var}'
    if cmax and (len(output) > cmax):
        output = output[:cmax] + ' ...'
    print(output)
    return output


def prn_list(aname, alist, m=3):
    outstring = f'> {aname} = [' + ', '.join(f'{a:.{m}f}' for a in alist) + ']; '
    # print(outstring)
    return outstring


def renormalize_x(aa, aminmax=(), bminmax=(), type=float):
    if len(aminmax) == 0:
        aminmax = (np.min(aa), np.max(aa))
    amin, amax = aminmax
    if len(bminmax) == 0:
        bminmax = aminmax
    bmin, bmax = bminmax
    aa_ = np.clip(aa, amin, amax)
    bb_ = (aa_ - amin) * (bmax - bmin) / (amax - amin) + bmin
    bb = bb_.astype(type)
    return bb


def find_param(text, varname, typ=int):
    """find numerical value of a variable: works for scalar or list of int or float"""

    text = str.lower(text)
    varname = str.lower(varname)
    t = text.find(varname)
    if t < 0:
        return np.nan

    txt = text[t:]
    txt = txt[txt.find('=')+1:]
    txt = txt[:txt.find(';')]
    # print(f'< txt = {txt}')

    t = txt.find('[')
    if t < 0:
        val = typ(float(txt))
    else:
        txt = txt[t+1: -1]
        # print(f'<1 txt = {txt}')
        vv = txt.split(',')
        val = [typ(float(v)) for v in vv]


    # aa = txt.split(',')
    # bb = [typ(float(a)) for a in aa]
    # if len(bb) == 1:
    #     bb = bb[0]

    # if typ == 'float': val = float(txt)
    # elif typ == 'int': val = int(txt)
    # else: val = txt
    # return val

    return val


# def find_param(text, varname, typ=int):
#     """find numerical value of a variable: works for scalar or list of int or float"""
#
#     text = str.lower(text)
#     varname = str.lower(varname)
#     t = text.find(varname)
#     if t < 0:
#         return np.nan
#
#     txt = text[t:]
#     txt = txt[txt.find('=')+1:]
#     txt = txt[:txt.find(';')]
#     # print(f'< txt = {txt}')
#
#     aa = txt.split(',')
#     bb = [typ(float(a)) for a in aa]
#     if len(bb) == 1:
#         bb = bb[0]
#
#     # if typ == 'float': val = float(txt)
#     # elif typ == 'int': val = int(txt)
#     # else: val = txt
#     # return val
#
#     return bb
#
#
# def find_paramlist(text, varname, typ=int):
#     """ for old data files """
#     text = str.lower(text)
#     varname = str.lower(varname)
#     t = text.find(varname)
#     if t < 0:
#         return np.nan
#
#     txt = text[t:]
#     txt = txt[txt.find('=')+1:]
#     txt = txt.replace('[', '').replace(']', ';')
#
#     aa = txt.split()


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
    from tkinter import filedialog

    # this_path = os.path.abspath(__file__)
    # print(runstamp(this_path))
    # f_path = filedialog.askopenfilename()
    # s_path, dict_parts = path_parts(f_path)
    # print(f'short_path = {s_path}')
    # print(f'fname = {dict_parts["fname"]}')

    text = 'lam_ns = [0.000, 0.602, 0.602, 0.601, 0.601, 0.600]; psi = 18.3456;'
    bb = find_param(text, 'lam_ns', float)
    b = find_param(text, 'psi')
    print(f'< bb = {bb}')
    print(f'< b = {b}')



