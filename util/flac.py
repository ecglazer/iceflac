#!/usr/bin/env python

import sys
try:
    import numpy as np
except ImportError:
    print('Error: Failed importing "numpy" module.')
    print('Please install the module and add it to PYTHONPATH environment variable.')
    sys.exit(1)

try:
    # python 2.7 or later
    from collections import Counter
except ImportError:
    try:
        # backport from http://code.activestate.com/recipes/576611
        # must be installed somewhere on your PYTHONPATH
        from counter import Counter
    except ImportError:
        pass


# the precision of saved file
doubleprecision = False
sizeofint = 4
if doubleprecision:
    sizeoffloat = 8
    default_dtype = np.double
else:
    sizeoffloat = 4
    default_dtype = np.single



class Flac(object):
    '''Read Flac data file. Most data are 2D arrays with shape (nx, nz) or (nex, nez).
    '''

    def __init__(self, swap_endian=False):
        self.swap_endian = swap_endian
        self.reload()
        return


    def reload(self):
        # read record number
        tmp = np.fromfile('../examples/_contents.0', sep=' ')
        tmp.shape = (-1, 3)
        self.time = tmp[:,2]
        self.nrec = len(self.time)

        # number of elements in x and z
        nex, nez = np.fromfile('../examples/nxnz.0', sep=' ', dtype=int)
        # number of nodes in x and z
        self.nx, self.nz = nex+1, nez+1
        self.nnodes = self.nx * self.nz
        self.nelements = nex * nez

        return


    def read_mesh(self, frame):
        columns = 2
        f = open('../examples/mesh.0')
        offset = (frame-1) * columns * self.nnodes * sizeoffloat
        f.seek(offset)
        x, z = self._read_data(f, columns)
        self._reshape_nodal_fields(x, z)
        return x, z


    def read_vel(self, frame):
        columns = 2
        f = open('../examples/vel.0')
        offset = (frame-1) * columns * self.nnodes * sizeoffloat
        f.seek(offset)
        vx, vz = self._read_data(f, columns)
        self._reshape_nodal_fields(vx, vz)
        return vx, vz


    def read_temperature(self, frame):
        columns = 1
        f = open('../examples/temperature.0')
        offset = (frame-1) * columns * self.nnodes * sizeoffloat
        f.seek(offset)
        T = self._read_data(f, columns)
        self._reshape_nodal_fields(T)
        return T


    def read_aps(self, frame):
        columns = 1
        f = open('../examples/aps.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        aps = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(aps)
        return aps


    def read_density(self, frame):
        columns = 1
        f = open('../examples/density.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        density = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(density)
        return density



    def read_eII(self, frame):
        columns = 1
        f = open('../examples/eII.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        eII = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(eII)
        return eII


    def read_sII(self, frame):
        columns = 1
        f = open('../examples/sII.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        sII = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(sII)
        return sII


    def read_sxx(self, frame):
        columns = 1
        f = open('../examples/sxx.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        sxx = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(sxx)
        return sxx


    def read_sxz(self, frame):
        columns = 1
        f = open('../examples/sxz.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        sxz = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(sxz)
        return sxz


    def read_szz(self, frame):
        columns = 1
        f = open('../examples/szz.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        szz = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(szz)
        return szz


    def read_srII(self, frame):
        columns = 1
        f = open('../examples/srII.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        srII = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(srII)
        return srII


    def read_pres(self, frame):
        columns = 1
        f = open('pres.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        pres = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(pres)
        return pres


    def read_diss(self, frame):
        columns = 1
        f = open('../examples/diss.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        diss = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(diss)
        return diss


    def read_visc(self, frame):
        columns = 1
        f = open('../examples/visc.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        visc = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(visc)
        return visc


    def read_phase(self, frame):
        columns = 1
        f = open('../examples/phase.0')
        offset = (frame-1) * columns * self.nelements * sizeofint
        f.seek(offset)
        phase = self._read_data(f, columns, count=self.nelements, dtype=np.int32)
        self._reshape_elemental_fields(phase)
        return phase


    def read_markers(self, frame):
        # read tracer size
        tmp = np.fromfile('_markers.0', sep=' ')
        tmp.shape = (-1, 4)
        n = int(tmp[frame-1,2])

        suffix = '.%06d.0' % frame
        dead = self._read_data('markdead' + suffix, count=n, dtype=np.int32)

        tmp = self._read_data('markx' + suffix, count=n)
        x = self._remove_dead_markers(tmp, dead)

        tmp = self._read_data('marky' + suffix, count=n)
        z = self._remove_dead_markers(tmp, dead)

        tmp = self._read_data('markage' + suffix, count=n)
        age = self._remove_dead_markers(tmp, dead)

        tmp = self._read_data('markphase' + suffix, count=n, dtype=np.int32)
        phase = self._remove_dead_markers(tmp, dead)

        tmp = np.arange(n)
        ID = self._remove_dead_markers(tmp, dead)

        return x, z, age, phase, ID


    def read_tracers(self):
        # read tracer size
        tmp = np.fromfile('_tracers.0', sep=' ')
        tmp.shape = (-1, 4)
        ntracerrec = tmp.shape[0]
        ntracers = int(tmp[0,1])
        n = ntracerrec * ntracers

        time = self._read_data('outtracktime.0', count=n)
        x = self._read_data('outtrackxx.0', count=n)
        x.shape = (ntracerrec, ntracers)

        z = self._read_data('outtrackyy.0', count=n)
        z.shape = (ntracerrec, ntracers)

        T = self._read_data('outtracktemp.0', count=n)
        T.shape = (ntracerrec, ntracers)

        p = self._read_data('outtrackpres.0', count=n)
        p.shape = (ntracerrec, ntracers)

        e = self._read_data('outtrackstrain.0', count=n)
        e.shape = (ntracerrec, ntracers)

        phase = self._read_data('outtrackphase.0', count=n)
        phase.shape = (ntracerrec, ntracers)

        return x, z, T, p, e, phase


    def _read_data(self, fileobj, columns=1,
                   count=None, dtype=None):
        '''Read data from a file-like object 'fileobj'.

        The 'dtype' specifies the storage type, default to single precision
        float.
        '''

        # number of nodes
        if count is None:
            count = self.nnodes

        # total number of items
        n = columns * count

        if dtype is None:
            dtype = default_dtype

        result = np.fromfile(fileobj, dtype, n)
        if self.swap_endian:
            result.byteswap(True)

        if columns == 1:
            return result
        else:
            result.shape = (columns, -1)
            return tuple(result[i,:] for i in range(columns))


    def _reshape_nodal_fields(self, *argv):
        for x in argv:
            x.shape = (self.nx, self.nz)
        return


    def _reshape_elemental_fields(self, *argv):
        for x in argv:
            x.shape = (self.nx-1, self.nz-1)
        return


    def _remove_dead_markers(self, a, dead):
        b = a[dead==1]
        return b


def elem_coord(x, z):
    '''Turning nodal coordinates to element coordinates'''
    cx = (x[:-1, :-1] + x[:-1, 1:] + x[1:, :-1] + x[1:, 1:]) / 4
    cz = (z[:-1, :-1] + z[:-1, 1:] + z[1:, :-1] + z[1:, 1:]) / 4
    return cx, cz


def make_uniform_grid(xmin, xmax, zmin, zmax, dx, dz):
    # grid size
    nx = (xmax - xmin) / dx + 1
    nz = (zmax - zmin) / dz + 1

    # generate uniform grid
    xx = np.linspace(xmin, xmax, nx)
    zz = np.linspace(zmin, zmax, nz)

    # the order of argument ensures the shape of arrays is (nx, nz)
    z, x = np.meshgrid(zz, xx)
    return x, z


def nearest_neighbor_interpolation2d(x0, z0, f0, x, z):
    '''Interpolating field f0, which is defined on (x0, z0)
    to a new grid (x, z) using nearest neighbor method'''

    if x0.shape != z0.shape:
        raise Exception('x0 and z0 arrays have different shape')

    if x0.shape != f0.shape:
        raise Exception('x0 and f0 arrays have different shape')

    if x.shape != z.shape:
        raise Exception('x and z arrays have different shape')

    # using 1d index for x0, z0, f0
    x0 = x0.flat
    z0 = z0.flat
    f0 = f0.flat

    dx = x[1,0] - x[0,0]
    dz = z[0,1] - z[0,0]

    nx, nz = x.shape
    f = np.zeros(x.shape)
    for i in range(nx):
        for j in range(nz):
            dist2 = ((x[i,j] - x0) / dx)**2 + ((z[i,j] - z0) / dz)**2
            ind = np.argmin(dist2)
            f[i,j] = f0[ind]

    return f


def gaussian_interpolation2d(x0, z0, f0, x, z):
    '''Interpolating field f0, which is defined on (x0, z0)
    to a new grid (x, z) using nearest neighbor method'''

    if x0.shape != z0.shape:
        raise Exception('x0 and z0 arrays have different shape')

    if x0.shape != f0.shape:
        raise Exception('x0 and f0 arrays have different shape')

    if x.shape != z.shape:
        raise Exception('x and z arrays have different shape')

    # using 1d index for x0, z0, f0
    x0 = x0.flat
    z0 = z0.flat
    f0 = f0.flat

    dx = 1.5 * (x[1,0] - x[0,0])
    dz = 1.5 * (z[0,1] - z[0,0])

    f = np.zeros(x.shape)
    g = np.zeros(x.shape)
    for i in range(len(x0)):
        weight = np.exp(-((x - x0[i]) / dx)**2 - ((z - z0[i]) / dz)**2)
        f += weight * f0[i]
        g += weight

    return f / g


def printing(*args, **kwd):
    '''printing(arg1, arg2, ..., stream=None)

    stream: None -- output to stdout
            filename -- output to file
            file-like object -- output to the object
    '''

    stream = kwd.get('stream')
    if stream == None:
        stream = sys.stdout
    elif isinstance(stream, (str, unicode)):
        filename = stream
        stream = open(filename, 'w')

    # flatten numpy arrays
    try:
        args = tuple(x.flat for x in args)
    except AttributeError:
        pass

    narg = len(args)
    fmt = '%.15e' + '\t%.15e'*(narg-1)

    for items in zip(*args):
        print >> stream, fmt % tuple(items)
    return



## This is an example on how to use this module.
## Before running this module, 'cd' to a directory containing the flac data.
if __name__ == '__main__':

    fl = Flac()

    # read the last record of the mesh and temperature
    x, z = fl.read_mesh(fl.nrec-1)
    T = fl.read_temperature(fl.nrec-1)

    # print (x, z, T) to screen
    printing(x, z, T)

    print('# time =', fl.time[fl.nrec-1], 'Myrs')
