import numpy as np
import matplotlib.pyplot as plt
import mayavi.mlab as ml

figx0 = 6.5
font0 = 'Consolas 10'


def plotAA(AA, figname='plotAA', caption='', ulimit=(), sxy=(1, 1), cmap='gray', pause=1):

    ny, nx = np.shape(AA)
    sx, sy = sxy
    figx = figx0 * sx
    figy = figx * (ny / nx) * sy
    figxy = (figx, figy)

    if ulimit:
        AA = np.clip(AA, ulimit[0], ulimit[1])

    plt.figure(figname, figsize=figxy, tight_layout=True)
    plt.clf()
    plt.imshow(AA, cmap=cmap, aspect='auto')
    plt.title(caption)

    plt_end(pause)


def graphB(B, figname='graphB', caption='', ulimit=(), xlimit=(), sxy=(1, 1), line='#1f77b4', pause=1):

    sy0 = 0.3

    sx, sy = sxy
    figx = figx0 * sx
    figy = figx * sy0 * sy
    figxy = (figx, figy)

    B = np.transpose(B)
    nx = len(B)

    plt.figure(figname, figsize=figxy, tight_layout=True)
    plt.clf()

    if xlimit == ():
        plt.plot(B, line)
    else:
        x1, x2 = xlimit
        xx = x1 + np.arange(nx) * (x2 - x1) / nx
        plt.plot(xx, B, line)

    plt.title(caption)
    plt.autoscale(enable=True, axis='x', tight=True)
    if ulimit:
        plt.ylim(ulimit)
    plt.grid(True)

    plt_end(pause)


def plotAAB(AA, figname='plotAAB', capA='', capB='', cmap='gray', line='#1f77b4',
            ulimit=(), roi=(), sxy=(1, 1), pause=1, crsr=True, ay=3):
    """
    ay = 3: vertical size of plotAA relative graphB
    """

    ''' plotAA '''
    ny, nx = np.shape(AA)
    sx, sy = sxy
    figx = figx0 * sx
    figya = figx * (ny / nx) * sy
    figyb = figya / ay
    figxy = (figx, figya + figyb)

    if ulimit:
        AA = np.clip(AA, ulimit[0], ulimit[1])

    plt.figure(figname, figsize=figxy, tight_layout=True)
    plt.clf()
    plt.subplot2grid((ay + 1, 1), (0, 0), rowspan=ay)
    plt.imshow(AA, cmap=cmap, aspect='auto')
    plt.title(capA)

    ''' graphB '''
    if len(roi) == 0:
        roi = (nx//2, ny//2, 10, 10)
    ix, iy, rx, ry = roi
    B = AA[iy, :]
    if crsr:
        put_cursor(roi)

    plt.subplot2grid((ay + 1, 1), (ay, 0), rowspan=ay)
    plt.plot(B, line)
    plt.title(capB)
    plt.autoscale(enable=True, axis='x', tight=True)
    if ulimit:
        plt.ylim(ulimit)
    plt.grid(True)

    plt_end(pause)


def put_cursor(roi):

    cline = 'yellow'
    croi = 'cyan'
    alpha = 0.75

    ix, iy, rx, ry = roi

    rx1 = ix - rx//2
    rx2 = ix + rx//2
    ry1 = iy - ry//2
    ry2 = iy + ry//2

    plt.axhline(y=iy, color=cline, alpha=alpha)
    # plt.axvline(x=ix, color=cline, alpha=alpha)
    plt.plot([rx1, rx2, rx2, rx1, rx1], [ry1, ry1, ry2, ry2, ry1], color=croi, alpha=alpha)


def plt_end(pause):

    if pause:
        plt.pause(pause)
    else:
        plt.show()


def mayaviAA(AA, figname='mayaviAA', caption='', view=(70, 20), ulimit=(), sxy=(1, 1), cmap='jet'):
    # import mayavi.mlab as ml

    figx0, figy0 = (700, 500)

    el, az = view
    sx, sy = sxy
    figxy = (figx0 * sx, figy0 * sy)

    if ulimit:
        AA = np.clip(AA, ulimit[0], ulimit[1])
    ml.figure(figname, size=figxy)
    ml.clf()
    ml.surf(AA, colormap=cmap, warp_scale='auto')
    ml.text(.02, .02, caption, width=len(caption)*.015)
    ml.view(elevation=el, azimuth=az)
    ml.show(stop=True)


def graph_many(graphs, figname='graph_many', col_row=(1, 1), sxy=(1, 1), line='#1f77b4', pause=1):

    sy0 = 0.25

    ncol, nrow = col_row
    sx, sy = sxy
    figx = ncol * figx0 * sx
    figy = nrow * figx0 * sy0 * sy
    figxy = (figx, figy)

    plt.figure(figname, figsize=figxy, tight_layout=True)
    plt.clf()

    for graph in graphs:
        uu, (irow, icol), caption, xlimit, ylimit = graph
        nx = len(uu)

        index = (irow - 0) * ncol + (icol + 1)
        plt.subplot(nrow, ncol, index)

        if xlimit == ():
            plt.plot(uu, line)
            plt.gca().xaxis.set_visible(False)
        else:
            x1, x2 = xlimit
            xx = x1 + np.arange(nx) * (x2 - x1) / nx
            plt.plot(xx, uu, line)
        plt.autoscale(enable=True, axis='x', tight=True)
        if ylimit:
            plt.ylim(ylimit)
        plt.title(caption)
        plt.grid(True)

    plt_end(pause)


if __name__ == '__main__':
    pi2 = 2 * np.pi
    nx, ny = (1920, 1080)
    xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
    kk = pi2 / 500
    roi = (300, 375, 100, 100)

    aa = 10 * np.sin(kk * xx) * np.sin(kk * yy) + 5 * (np.random.rand(ny, nx) - .5)
    b = aa[roi[1], :]

    graphs = []
    graphs += [(aa[0, :], (0, 0), 'aa[0, :], (0, 0)', (), ())]
    graphs += [(aa[1, :], (1, 0), 'aa[1, :], (1, 0)', (), ())]
    graphs += [(aa[2, :], (2, 0), 'aa[2, :], (2, 0)', (), ())]
    graphs += [(aa[3, :], (1, 1), 'aa[3, :], (1, 1)', (), ())]
    graphs += [(aa[4, :], (2, 2), 'aa[4, :], (2, 2)', (), ())]
    graphs += [(aa[5, :], (4, 2), 'aa[5, :], (4, 2)', (0, nx), ())]

    plotAA(aa, pause=1)
    graphB(b, pause=1)
    plotAAB(aa, roi=roi, pause=1)
    graph_many(graphs, col_row=(3, 5), sxy=(.75, .75), pause=1)
    mayaviAA(aa)

    plt.show()


