# code adapted from Mingzhou Chen's Matlab code
import sys
from typing import List, Tuple

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from tqdm import trange


def cart2pol(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Summary

    Parameters
    ----------
    x : np.ndarray
        Description
    y : np.ndarray
        Description

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Description
    """

    theta = np.arctan2(y, x)
    rho = np.sqrt(x**2 + y**2)

    return theta, rho


def create_bottle(x: np.ndarray, y: np.ndarray,
                  centre: List[float], radius: float,
                  width: float, dx: float,
                  dz: float) -> np.ndarray:
    """set air to 2, set glass to 1, set alcohol to 0

    Parameters
    ----------
    x : np.ndarray
        Description
    y : np.ndarray
        Description
    centre : List[float]
        Description
    radius : float
        Description
    width : float
        Description
    dx : float
        Description
    dz : float
        Description

    Returns
    -------
    np.ndarray
        Description

    """

    real_width = width / dx
    real_radius = radius / dx
    real_centre = [centre[0]*dx, centre[1]*dx]

    out = np.ones_like(x)
    rin = np.sqrt((x - real_centre[0])**2 + (y - real_centre[1])**2)
    out[rin >= real_radius] = 2
    out[rin <= real_radius-real_width] = 0
    return out


def create_bottle2(x, y, centre, radius, width, dx, dz):

    out = np.ones((512, 512))
    real_centre = [centre[0]*dx, centre[1]*dz]

    x = 0
    for i in range(512):
        z = 0
        for j in range(512):
            if np.sqrt((z - real_centre[1])**2+(x - real_centre[0])**2) >= radius:
                out[i, j] = 2.
            elif np.sqrt((z - real_centre[1])**2+(x - real_centre[0])**2) <= radius - width:
                out[i, j] = 0.
            z += dz
        x += dx

    return out


def filled_space_prop(e: np.ndarray, arg: np.ndarray) -> np.ndarray:
    """Propagate beam through filled space


    Parameters
    ----------
    e : np.ndarray
        Electirc field
    arg : np.ndarray
        Distance to propagate
    k1k22k : np.ndarray
        k-space constant

    Returns
    -------
    out : np.ndarray
        New electric field
    """

    freq = np.exp(1j*arg)
    out = np.fft.ifft2(np.fft.fft2(e)*freq)

    return out


def free_space_prop(e: np.ndarray, distance: float, k1k22k: np.ndarray) -> np.ndarray:
    """Propagate beam through free space


    Parameters
    ----------
    e : np.ndarray
        Electirc field
    distance : float
        Distance to propagate
    k1k22k : np.ndarray
        k-space constant

    Returns
    -------
    out : np.ndarray
        New electric field
    """

    arg = distance*k1k22k
    freq = np.exp(1j*arg)
    out = np.fft.ifft2(np.fft.fft2(e)*freq)

    return out


w0 = 582*4  # um
wavelength = .785  # um
k = 2 * np.pi / wavelength
axicon_angle = 5  # degrees
n = 1.45  # refractive index of bottle

k_r = k * (n - 1) * (axicon_angle) * np.pi / 360.
ell = 0  # what is this?

xymax = 5000  # um
nxy = 512  # number of voxel in xy planes
nz = 1000  # number of voxels in z plane
zmax = w0 * (k / k_r)
L = 3 * zmax
R = L

dz = L / nz
xmax = xymax
dx = xmax / nxy
kmax = 2*np.pi/dx
dk = kmax / nxy
nmid = int(nxy / 2)

v = np.arange(0, nxy)
x, y = np.meshgrid(v, v)
x = x*dx - xmax/2
y = y*dx - xmax/2

p = v[v > nmid]
v[p] = nxy - v[p]
v = v * dk

k2, k1 = np.meshgrid(v, v)
k1k22k = -dz*(k1**2 + k2**2)/(2.*k)
theta, r = cart2pol(x, y)

e = np.exp(-(r-1612)**2/300**2)  # magic numbers? possible radius and width?
e = np.abs(e) * np.exp(1j*ell*theta)
iprof = []
lenses = []
ztotal = 0
print("Free Space before lens")
for i in trange(int(nz/10)):
    e = free_space_prop(e, 1., k1k22k)
    absy = np.abs(e[nmid, :])**2
    imax = np.amax(absy)
    iprof.append(absy / imax)
    lenses.append(np.zeros_like(absy))
    # ztotal += dz
# e = free_space_prop(e, dz*nz, k1k22k)

# propagate to bottles edge after focusing lens
e *= np.exp(-1j*k*r**2/(2.*R))  # lens
# propagate to bottles surface
print("Free Space after lens before bottle")
for i in trange(int(nz/2)+0):
    # e = free_space_prop(e, 1., k1k22k)
    absy = np.abs(e[nmid, :])**2
    imax = np.amax(absy)
    iprof.append(absy / imax)
    lenses.append(np.zeros_like(absy))
    ztotal += dz


# print(160061, 40000 / dz, R)
# sys.exit(0)
dzz = 40./nz
dxx = 3./300
L0 = 2  # bottlethickness, mm
RR = 40  # bottle radius, mm
d = 0
bottlePS = np.exp(-1j*k*1.5*r**2/2/(RR/40*R))
bottlePS = np.stack([bottlePS[:, int(nxy/2)]for j in range(512)])
# e *= bottlePS

print("In bottle")
for i in trange(385):
    # roi = np.zeros_like(k1k22k)
    # d += dzz/2
    # L1 = 0
    # if d <= L0:
    #     LL = np.sqrt(RR**2 - (RR-d)**2)
    #     LL = int(round(LL / dxx))
    # else:
    #     L1 = np.sqrt((RR-L0)**2-(RR-d)**2)
    #     LL = np.sqrt(RR**2-(RR-d)**2)
    #     L1 = int(round(L1/dxx))
    #     LL = int(round(LL/dxx))

    # d = d + dzz/2
    # if L1 == 0:
    #     if LL >= nxy/2:
    #         roi = roi + 1
    #     else:
    #         roi[int(nxy/2-LL):int(nxy/2+LL), 1:nxy] = 1
    # else:
    #     if LL < nxy/2:
    #         roi[int(nxy/2)+L1:int(nxy/2)+LL, 1:nxy] = 1
    #         roi[int(nxy/2)-LL:int(nxy/2)-L1, 1:nxy] = 1
    #     elif LL > nxy/2 and L1 < nxy/2:
    #         roi[int(nxy/2)+L1:nxy, 1:nxy] = 1
    #         roi[1:int(nxy/2)-L1, 1:nxy] = 1
    # roi = np.fft.fftshift(roi)
    # argg = -dz*(k1**2+k2**2)/2/k/n
    # if d >= L0:
    #     arg1 = -dz*(k1**2+k2**2)/2/k/1.3
    # else:
    #     arg1 = -dz*(k1**2+k2**2)/2/k

    # argnew = argg*roi+arg1*(1-roi)
    # freqnew = np.exp(1j*argnew)  # new propagator

    # e = np.fft.fft2(np.fft.ifft2(e)*freqnew)
    absy = np.abs(e[nmid, :])**2
    imax = np.amax(absy)
    iprof.append(absy / imax)
    lenses.append(np.zeros_like(absy))
    ztotal += dz

# save file
out_int = np.abs(e.T)**2
out_int.tofile("bessel-normal.dat")
print(ztotal, len(iprof), dz, len(iprof)*dz)
fig, axs = plt.subplots(1, 2)
extent = [0, ztotal*1e-3, 1e-3*-512*dx/2., 1e-3*512*dx/2.]
# extent = [0-(1e-3)*(nz*dz)/40, (1e-3*dz*len(iprof)/4)-(1e-3*nz*dz)/10, 1e-3*-512*dx/2., 1e-3*512*dx/2.]
axs[0].imshow(np.array(iprof).T, extent=extent, interpolation="gaussian", aspect="auto")
axs[0].set_xlabel("z/mm")
axs[0].set_ylabel("x/mm")

axs[1].imshow(np.abs(e[115:395, 115:395].T)**2, interpolation="gaussian")
plt.show()
