#!/usr/bin/env python
'''
Read flac density and coordinate data and return gravity and topography at
regular spacing.
Filling basin with sediment (to 2km below sea level) and erode forearc height
to sea level before computing gravity.
'''

import sys
import numpy as np
import flac

gn = 6.6726e-11
ocean_density = 1030.
sediment_density = 2400.

def find_trench_index(zz):
    '''Returns the i index of trench location.'''
    # the highest point defines the forearc
    imax = zz.argmax()
    # the trench is the lowest point west of forearc
    i = zz[:imax].argmin()
    return i


def find_peaks(zz):
    nx = zz.shape[0]

    # max basin depth
    basin_min_depth = -2000
    itrench = find_trench_index(zz)

    peaks = []
    for i in range(itrench+1, nx-1):
        left_dz = zz[i] - zz[i-1]
        right_dz = zz[i+1] - zz[i]
        if left_dz * right_dz <= 0 and (left_dz > 0 or right_dz < 0):
                peaks.append(i)

    # add the right boundary as a peak
    peaks.append(nx-1)
    return peaks



def compute_gravity(frame):
    ## read data
    fl = flac.Flac()
    x, z = fl.read_mesh(frame)  # in km
    x *= 1e3
    z *= 1e3

    # surface coordinates
    xx = x[:,0]
    zz = z[:,0]
    xmin = x[0,0]
    xmax = x[-1,0]

    # center of elements
    cx, cz = flac.elem_coord(x, z)

    # area of elements = 0.5 * | AC x BD | = 0.5 * |xdiag1*zdiag2 - xdiag2*zdiag1|
    #  A -- D
    #  |    |
    #  B -- C
    xdiag1 = x[0:-1, 0:-1] - x[1:, 1:]
    zdiag1 = z[0:-1, 0:-1] - z[1:, 1:]
    xdiag2 = x[1:, 0:-1] - x[0:-1, 1:]
    zdiag2 = z[1:, 0:-1] - z[0:-1, 1:]
    area = 0.5 * np.abs(xdiag1*zdiag2 - xdiag2*zdiag1)


    rho = fl.read_density(frame)  # in kg/m^3

    # anything above sea level is removed
    rho[cz > 0] = 0

    ## benchmark case: infinite-long cylinder with radius R
    ## buried at depth D
    #R = 10e3
    #D = -150e3
    #drho = 1000
    #rho = np.zeros(cx.shape)
    #midx = 0.5 * (xmin + xmax)
    #midz = 0.5 * z.min()
    #dist2 = (x - midx)**2 + (z - D)**2
    #rho[dist2 < R**2] = drho
    #ocean_density = 0

    mass = rho * area


    ## calculate gravity at these points
    # px in uniform spacing
    px = np.linspace(xmin, xmax, num=5*fl.nx)
    # pz is a few km above the highest topography to avoid high frequency oscillation
    pz_height = max(0, np.max(zz)) + 4e3
    print 'gravity evaluated at %f km' % pz_height
    pz = np.ones(px.shape) * pz_height

    # original topography defined in px grid
    topo = np.interp(px, x[:,0], z[:,0])

    ## contribution of material inside the model
    grav = np.empty(px.shape)
    for i in range(len(grav)):
        dx = px[i] - cx
        dz = pz[i] - cz
        dist2 = (dx**2 + dz**2)
        # downward component of gravity of an infinite long line source of density anomaly
        # see Turcotte & Schubert, 1st Ed., Eq 5-104
        grav[i] = 2 * gn * np.sum(mass * dz / dist2)


    ## contribution of sedimentary basin, only to the right of trench
    peaks = find_peaks(zz)
    itrench = find_trench_index(zz)
    basin_depth = -2000
    sed_density = 2200
    sed_thickness = np.zeros(fl.nx)

    for ii in range(len(peaks)-1):
        fill_height = min((basin_depth, zz[peaks[ii]], zz[peaks[ii+1]]))
        for i in range(peaks[ii], peaks[ii+1]):
            if zz[i] < fill_height:
                sed_thickness[i] = fill_height - zz[i]
                zz[i] = fill_height

    for i in range(itrench, fl.nx-1):
        sedz = 0.5 * (sed_thickness[i] + sed_thickness[i+1])
        if sedz > 0:
            midz = 0.5 * (zz[i] + zz[i+1])
            midx = 0.5 * (xx[i] + xx[i+1])
            m = (xx[i+1] - xx[i]) * sedz * sed_density
            dx = px - midx
            dz = pz - midz
            dist2 = (dx**2 + dz**2)
            grav += 2 * gn * m * dz / dist2


    ## contribution of ocean
    for i in range(fl.nx-1):
        midz = 0.5 * (zz[i] + zz[i+1])
        if midz < 0:
            midx = 0.5 * (xx[i] + xx[i+1])
            m = (xx[i+1] - xx[i]) * -midz * ocean_density
            dx = px - midx
            dz = pz - midz
            dist2 = (dx**2 + dz**2)
            grav += 2 * gn * m * dz / dist2

    # contribution of material outside left boundary
    # assuming the leftmost element extend to negative infinity
    # see Turcotte & Schubert, 1st Ed., Eq 5-106
    for i in range(fl.nz-1):
        sigma = rho[0,i] * (z[0,i] - z[0,i+1])
        dx = px - xmin
        dz = pz - 0.5 * (z[0,i] + z[0,i+1])
        angle = np.arctan2(dx, dz)
        grav += 2 * gn * sigma * (0.5*np.pi - angle)

    if zz[0] < 0:
        sigma = ocean_density * -zz[0]
        dx = px - xmin
        dz = pz - 0.5 * zz[0]
        angle = np.arctan2(dx, dz)
        grav += 2 * gn * sigma * (0.5*np.pi - angle)


    # ditto for right boundary
    for i in range(fl.nz-1):
        sigma = rho[-1,i] * (z[-1,i] - z[-1,i+1])
        dx = xmax - px
        dz = pz - 0.5 * (z[-1,i] + z[-1,i+1])
        angle = np.arctan2(dx, dz)
        grav += 2 * gn * sigma * (0.5*np.pi - angle)

    if zz[-1] < 0:
        sigma = ocean_density * -zz[-1]
        dx = xmax - px
        dz = pz - 0.5 * zz[-1]
        angle = np.arctan2(dx, dz)
        grav += 2 * gn * sigma * (0.5*np.pi - angle)

    ## set reference gravity
    ## far right gravity to 0
    grav -= grav[-1]

    return px, topo, grav


def modify_topo(topo):
    peaks = find_peaks(topo)
    itrench = find_trench_index(topo)
    basin_depth = -2000

    topomod = topo.copy()

    # erode forearc height to sea level
    topomod[topo > 0] = 0

    # fill basin
    for ii in range(len(peaks)-1):
        fill_height = min((basin_depth, topo[peaks[ii]], topo[peaks[ii+1]]))
        for i in range(peaks[ii], peaks[ii+1]):
            if topo[i] < fill_height:
                topomod[i] = fill_height

    return topomod


def compute_gravity2(frame):
    '''Similar to compute_gravity(), but also return modified topo'''
    px, topo, grav = compute_gravity(frame)
    topomod = modify_topo(topo)

    return px, topo, topomod, grav


if __name__ == '__main__':
    frame = int(sys.argv[1])

    px, topo, topomod, gravity = compute_gravity2(frame)
    flac.printing(px, topo, topomod, gravity)


