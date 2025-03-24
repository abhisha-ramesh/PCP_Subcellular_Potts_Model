# A Subcellular Potts Model for Planar Cell Polarity 
This project simulates Planar Cell Polarity using Cellular Potts Model (aka Glazier-Graner-Hogeweg Model) using the CompuCell3D(CC3D) (Version 4.2.3) simulation software. The model investigates conditions under which tissue-wide polarity can be established in the absence of a long-range morphogen gradient. This simulation is part of the research presented in the paper "Start Small: A Sub-cellular Potts Model for Tissue-wide Planar Cell Polarity without Morphogens" by Abhisha Thayambath and Julio M. Belmonte. This README file provides instructions on how to setup, run and modify this simulation using CC3D simulation environment.

## Installing CC3D
The latest version of CC3D can be downloaded and installed from the [official website](https://compucell3d.org/) (supported on Windows, Mac and Linux). Simulations can be executed using [CompuCell3D Player](https://github.com/CompuCell3D/cc3d-player5/tree/master).

## How to run and modify simulation?
The folder PCP_code contains a .cc3d project file and a directory called Simulation where one can find all the Python Scripts. To load a simulation in CC3D Player or into [Twedit++](https://github.com/CompuCell3D/cc3d-twedit5/tree/master) (IDE for CC3D), use the .cc3d file. The simulation can be modified by loading the python scripts into Twedit++.
