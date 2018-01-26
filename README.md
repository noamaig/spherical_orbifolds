# Spherical Orbifold Tutte Embeddings

Matlab code implementing the [Siggraph 2017 paper, "Spherical Orbifold Tutte Embeddings"](https://noamaig.github.io/html/projects/spherical/spherical_low.pdf).

An extension of Tutte's embedding to spherical domains, namely spherical orbifolds - tilings of the sphere.


### Example scripts:
- `script_orbifold_sphere` - map a sphere-topology mesh to a spherical orbifold.
- `script_orbifold_disk` - map a disk-topology mesh to a spherical orbifold.
- `script_NOT_ORBIFOLD_spherical_tutte` - embed a disk-topology mesh into a convex spherical polygon.

### Main files to look at:
- `embed_from_data` - compute the spherical embedding from a given mesh and choice of cone points. The detection of sphere\disk topology is automatic. If you want to embed into a convex polygon the input should be the positions of the boundary instead of the cone data.
- `Embedding` - the object representing the computed embedding. Has all the data as well as functions for visualization (see the drawing part in the example scripts).

The code is provided as-is for academic use only and without any guarantees. Please contact the author to report any bugs.
Written by [Noam Aigerman](https://noamaig.github.io/) and [Shahar Kovalsky](https://github.com/shaharkov).

 
