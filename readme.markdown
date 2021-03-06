Summer 2011
===========

A collection of code developed over summer 2011 primarily by Allan Liu.
This code will be integrated in time into other packages/projects.


Elliptic Integrals
------------------

Igor Moiseev has a collection of Matlab code for calculating
elliptic integrals and assorted special functions.
The official code repository is on Google Code: <http://code.google.com/p/elliptic/>

We will be extending the functions `elliptic12` and `elliptic3` to accept
greater ranges of parameter input as well as computing the complete elliptic
integrals for these greater ranges.

Files:

* `elliptic123.m`
* `elliptic_tests.m`


Magnet code
-----------

We will be implementing extensions to the Matlab code for calculating forces,
etc., between magnets that is located on GitHub: <http://github.com/wspr/magcode>

Files:

* `Torque.m` — Calculating the torque between two parallel cuboid magnets.
* `magnetcoil.m` — Calculating the axial between a coaxial magnet and solenoid.

ANSYS
-----

Various scripts for calculating force and torque for a variety of geometries
in Ansys. See the files in the relevant folder for further information.

Copyright and licensing
-----------------------

The code in this repository is distributed under the Apache License (v2)
for now; its licensing may change when contributed into new projects. The
Apache License is located at: <http://www.apache.org/licenses/LICENSE-2.0>

Copyright 2011 Allan Liu and Will Robertson
