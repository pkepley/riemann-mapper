import os, sys

sys.path.append(os.getcwd() + "/../")

import numpy as np
from svgpathtools import svg2paths, disvg, path
from PointsOnGamma import PointsOnGamma
from KernelsKS import KernelsKS
from RiemannMapper import RiemannMapper
import matplotlib.pyplot as plt

save_image_flag = True

# Indiana image
paths, attrs = svg2paths("../imgs/indiana_map.svg")

# Grab the path. Note that it's upside down, so we'll correct that first
gamma = paths[0].rotated(180)

# Distribute n points around a piece-wise continuous path,
n = 5000
pog = PointsOnGamma(gamma, n)
zs_b = pog.ps
xs_b = zs_b.real
ys_b = zs_b.imag

# Get a grid of points to plot
xlo, xhi = np.min(xs_b), np.max(xs_b)
ylo, yhi = np.min(ys_b), np.max(ys_b)
nx, ny = 100, 100
nps = nx * ny
xs, ys = np.linspace(xlo, xhi, nx), np.linspace(ylo, yhi, ny)
ps = np.array(np.meshgrid(xs, ys))
xs, ys = ps[0, :, :].reshape((nps, 1)), ps[1, :, :].reshape((nps, 1))
ps = xs + (1.0j) * ys

# Determine if a point is inside of Omega by checking if the
# winding number is 1
winding_at_ps = pog.cauchy_integral(np.ones(zs_b.shape), ps.T).reshape((nps, 1))
print(np.abs(winding_at_ps))
i_xs = np.abs(winding_at_ps - 1) < 0.2
ws = ps[i_xs]

# Plot the path and grid
plt.figure(1, figsize=(12, 12))
ws_xx, ws_yy = ws.real, ws.imag
plt.plot(xs_b, ys_b)
plt.scatter(ws_xx, ws_yy, marker=".")
plt.axes().set_aspect("equal")

# not guaranteed to fall inside of gamma, but set a the mean of gamma
c = np.mean(ps)
v = ps[n // 4] - c
a = c + 0.0 * v

# Create the Kurzmann-Stein
kks = KernelsKS(a, pog)
kks.solve_ks()

rm = RiemannMapper(kks)
zz = rm.evaluate_riemann(ws)

x_zz, y_zz = zz[0].real, zz[0].imag
plt.figure(2, figsize=(12, 12))
tt = np.linspace(0, 2 * np.pi, 100)
circ_xs, circ_ys = np.cos(tt), np.sin(tt)
plt.plot(circ_xs, circ_ys, label="circle")
plt.scatter(x_zz, y_zz, marker=".")
plt.axis([-1, 1, -1, 1])
plt.axes().set_aspect("equal")

if save_image_flag:
    plt.figure(1)
    plt.savefig("../imgs/indianaSource.png", bbox_inches="tight")
    plt.figure(2)
    plt.savefig("../imgs/indianaTarget.png", bbox_inches="tight")

plt.show()
