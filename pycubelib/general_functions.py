import numpy as np
import datetime
import os

pi = np.pi
pi2 = pi * 2


def runstamp(script_path):
    s_path, _ = path_parts(script_path, 3)
    time_stamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    stamp = f'[{s_path}] {time_stamp}'
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


def find_param(text, varname, typ=int):
    """find numerical value of a variable and convert to typ"""

    text = str.lower(text)
    varname = str.lower(varname)
    t1 = text.find(varname)
    if t1 < 0:
        return np.nan

    txt = text[t1:]
    txt = txt[txt.find('=')+1:]
    txt = txt[:txt.find(';')]

    t1 = txt.find('[')
    if t1 < 0:
        val = typ(float(txt))
    else:
        t2 = txt.find(']')
        txt = txt[t1+1: t2-1]
        vv = txt.split(',')
        val = [typ(float(v)) for v in vv]

    return val


if __name__ == '__main__':

    text = 'lam_ns = [0.000, 0.602, 0.602, 0.601, 0.601, 0.600]; psi = 18.3456;'
    bb = find_param(text, 'lam_ns', float)
    b = find_param(text, 'psi')
    print(f'< bb = {bb}')
    print(f'< b = {b}')



