This folder is for tools created for modeling purposes.

1) checkneighbors.cpp

This file will read a lammps file and find neighbors around every atom. Uses boost library and is very fast.
This is useful to create outer shell of structures from solid chunks or blocks or complex geometries.

compile using
 g++ -o executable.out checkneighbors.cpp -lboost_system
