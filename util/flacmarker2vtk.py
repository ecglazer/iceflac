#!/usr/bin/env python

'''Convert the binary marker output of flac to VTK (vtp) files.
'''

import sys, os
import numpy as np
import flac
from flac2vtk import vts_dataarray


# filtering markers to only those within the domain bounds (in km)
filtering = True
xmin = 300
xmax = 700
zmin = -50
zmax = 100


def filter_marker(x, z, age, phase, ID):
    # bool * bool is element-wise logical AND
    ind = (xmin <= x) * (x <= xmax) * (zmin <= z) * (z <= zmax)
    x = x[ind]
    z = z[ind]
    age = age[ind]
    phase = phase[ind]
    ID = ID[ind]
    return x, z, age, phase, ID


def main(path, start=1, end=-1):

    # changing directory
    os.chdir(path)

    fl = flac.Flac()
    if end == -1:
        end = fl.nrec

    for i in range(start, end+1):
        x, z, age, phase, ID = fl.read_markers(i)
        if filtering:
            x, z, age, phase, ID = filter_marker(x, z, age, phase, ID)
        nmarkers = len(x)

        print 'Writing record #%d, model time=%.3e, %d markers' % (i, fl.time[i-1], nmarkers)
        fvtp = open('flacmarker.%06d.vtp' % i, 'w')
        vtp_header(fvtp, nmarkers)

        # point-based data
        fvtp.write('  <PointData>\n')
        vts_dataarray(fvtp, age, 'age', 1)
        vts_dataarray(fvtp, phase.astype(np.int32), 'phase', 1)
        vts_dataarray(fvtp, ID.astype(np.int32), 'ID', 1)
        fvtp.write('  </PointData>\n')

        # point coordinates

        # VTK requires vector field (velocity, coordinate) has 3 components.
        # Allocating a 3-vector tmp array for VTK data output.
        tmp = np.zeros((nmarkers, 3), dtype=x.dtype)
        tmp[:,0] = x
        tmp[:,1] = z
        fvtp.write('  <Points>\n')
        vts_dataarray(fvtp, tmp, '', 3)
        fvtp.write('  </Points>\n')

        vtp_footer(fvtp)
        fvtp.close()
    return


def vtp_header(f, npoints):
    f.write(
'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
<PolyData>
<Piece NumberOfPoints="{0}">
'''.format(npoints))
    return


def vtp_footer(f):
    f.write(
'''</Piece>
</PolyData>
</VTKFile>
''')
    return


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print '''usage: flacmarker2vtk.py path [step_min [step_max]]

Processing flac marker output to VTK format.
If step_max is not given, processing to latest steps
If both step_min and step_max are not given, processing all steps'''
        sys.exit(1)

    path = sys.argv[1]

    start = 1
    end = -1
    if len(sys.argv) >= 3:
        start = int(sys.argv[2])
        if len(sys.argv) >= 4:
            end = int(sys.argv[3])

    main(path, start, end)
