# Simple 3D Computational Geometry

SCG is a simple, value oriented, MIT License library for 3D computational geometry in C++11.

## Examples

## Finding the angle between two vectors

```cpp
point p(1, 1, 0);
point q(-1, -1, 0);
double a = angle_between(p, q);
assert(within_eps(a, 180, 0.00001));

```

### Finding a rotation from point vector onto another

```cpp
point from(1, 0, 0);
point to(0, 0, 1);

const rotation r = rotate_from_to(from, to);

cout << apply(r, from) << endl;
```

### Loading a mesh from an STL file

```cpp
auto triangles = parse_stl("./test/stl-files/SlicedCone.stl").triangles;
triangular_mesh m = make_mesh(triangles, 0.001);
cout << "# of triangles = " << m.triangle_list().size() << endl;
```
## Build Instructions

After cloning the repo cd into the top directory of the project. Then do:

```bash
cmake .
make -j
```

This will build several things.
* libutils.a: A utility library used by the main geometry library
* libgeometry.a: The geometry library itself.
* geometry-tests: The unit tests for this library
* run-examples: An executable that contains some samples of how to use this library

To run the tests do:

```bash
./geometry-tests 
```

To run the examples do:

```bash
./run-examples
```

## Data Structures

* 3D points
* Planes
* Triangles
* Boxes
* Triangular meshes (stored using a halfedge data structure)
* Rotations and homogeneous transformations
* Height fields
* Voxel volumes

SCG also has a built in binary STL file parser to load data in to triangular meshes.

## Dependecies
* boost
