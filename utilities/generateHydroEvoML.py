#!/usr/bin/env python3
# coding: utf-8

import numpy as np


def outputEvoInMUSICFormat(evoData):
    filename = "evolution_all_xyeta.dat"
    f = open(filename, 'wb')

    nchannels, nx, ny, ntau = evoData.shape
    tau0 = 0.4
    dtau = 0.1
    dx = 0.265625
    dy = 0.265625
    xmin = -17.
    ymin = -17.
    neta = 1
    deta = 0.1
    etamin = 0.0
    header = [tau0, dtau, nx, dx, xmin, ny, dy, ymin, neta, deta, etamin,
              0., 1., 1., 0., 11]
    f.write(np.array(header, dtype='float32').tobytes())

    ieta = 0
    eta_s = 0.
    cosh_eta = np.cosh(eta_s)
    sinh_eta = np.sinh(eta_s)
    cs2 = 1./3.
    for itau in range(ntau):
        for ix in range(nx):
            for iy in range(ny):
                T = evoData[0, ix, iy, itau]
                e = evoData[1, ix, iy, itau]
                P = evoData[2, ix, iy, itau]
                ux = evoData[3, ix, iy, itau]
                uy = evoData[4, ix, iy, itau]
                ueta = evoData[5, ix, iy, itau]
                utau = np.sqrt(1. + ux**2 + uy**2 + ueta**2)
                uz = ueta*cosh_eta + utau*sinh_eta
                ideal = [itau, ix, iy, ieta, e, P, T, cs2, ux, uy, uz]
                f.write(np.array(ideal, dtype='float32').tobytes())
    f.close()


data_modeled = np.load(f"FNO-outputs/data_modeled.npy")

eventId = 23
outputEvoInMUSICFormat(data_modeled[eventId])


