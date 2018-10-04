# Simple 3D Computational Geometry

SCG is a simple, value oriented, MIT License library for 3D computational geometry in C++11.

## Data Structures

* 3D points
* Planes
* Triangles
* Triangular meshes (stored using a halfedge data structure)
* Rotations and homogeneous transformations
* Height fields
* Voxel volumes

SCG also has a built in binary STL file parser to load data in to triangular meshes.

## Examples

### Computing the rotation that will map one 3D vector onto another

```cpp
point from(1, 1, 0);
point to(0, 0, 1);

const rotation r = rotate_from_to(from, to);

```

### Loading a mesh from an STL file

```cpp
auto triangles = parse_stl("./test/stl-files/SlicedCone.stl").triangles;
triangular_mesh m = make_mesh(triangles, 0.001);
```

## Dependecies
* boost
