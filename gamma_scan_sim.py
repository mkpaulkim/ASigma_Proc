import numpy as np
import matplotlib.pyplot as plt
pi = np.pi
pi2 = pi * 2

ax = 10000
nx = 2000
hh = 2000
lam1 = 0.60180057
lam2 = 0.60172814
nw = 20
z0 = 50

dlam = lam2 - lam1
zz = np.arange(nx) * hh / nx - hh/2

lams = [0]
lams1n = [0]
eee = np.zeros((nw, nx)) * 0j
for n in range(1, nw+1):
    lam = lam1 + (n-1) * dlam
    lams += [lam]
    if n == 1:
        lam1n = 0
    else:
        lam1n = lam1 * lam / (lam1 - lam)
    print(f'> n = {n}: lam = {lam:.8f}; lam1n = {lam1n:.1f}')
    lams1n += [lam1n]

    kn = pi2 / lam
    ee = np.exp(1j * kn * zz) * np.exp(1j * kn * z0)
    eee[n-1, :] = ee

fff = np.fft.fft(eee, axis=0)
fff = np.fft.fftshift(fff, axes=0)
ffa = np.abs(fff)

lam12 = lams[1] * lams[2] / (lams[1] - lams[2])
z_step = lam12 / (nw - 1)
zz_ffa = np.argmax(ffa, axis=0) * z_step - lam12/2

plt.figure(1, tight_layout=True)
plt.plot(zz)
plt.autoscale(enable=True, axis='x', tight=True)

plt.figure(2, tight_layout=True)
plt.imshow(ffa, aspect='auto', origin='lower', cmap='gray')

plt.figure(3, tight_layout=True)
plt.plot(zz_ffa)
plt.autoscale(enable=True, axis='x', tight=True)

plt.show()






