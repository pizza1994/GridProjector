# GridProjector

This code removes the external polyhedra of a conformal grid and projects the result on the target surface. 
## Dependencies

This code is built on top of [Cinolib](https://github.com/mlivesu/cinolib.git).

## Building
Clone this repository, including submodules, with:
```
git clone --recursive https://github.com/pizza1994/GridProjector.git
```
Build the executable by running the following commands:
```
cd GridProjector 
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=<build type> ..
make
```
## Usage
You can project a conformal grid on a target surface by running the command below. Substitute ``<projection_method>`` with mesh_smoother or grid_projector
```
./GridProjector <surface_mesh.obj> <conformal_grid.mesh> <output.mesh> <projection_method>
```