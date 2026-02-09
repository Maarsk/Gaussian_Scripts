
import ase.io.cube

import numpy as np

import os



fname = 'MOvalue.cub'
print("Reading cube file {}".format(fname))
data, atoms = ase.io.cube.read_cube_data(fname)

# Here, I want the electron density, not the wave function
data = data**2

# If data is too large, just reduce it by striding with steps >1
sx, sy, sz = 1, 1, 1
data = data[::sx,::sy,::sz]

# Note the reversed order!!
nz, ny, nx = data.shape
nframes = 1
header = np.array([nx,ny,nz,nframes])

#open and write to file
vdata = data.flatten() / np.max(data)
vfname = os.path.splitext(fname)[0] + '.bvox'
vfile = open(vfname,'wb')
print("Writing Blender voxel file {}".format(vfname))
header.astype('<i4').tofile(vfile)
vdata.astype('<f4').tofile(vfile)