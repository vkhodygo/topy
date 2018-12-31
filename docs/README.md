The docs folder contains [William Hunter master's dissertation](thesis.pdf) since
ToPy was part of it. The maths behind topology optimization is explained in it.
There is also quite a bit about ToPy in it, but now a little outdated.
The persistent URL for the dissertation is: http://hdl.handle.net/10019.1/2648

Other useful documentation pertaining to ToPy and/or topology optimisation
has been added to this folder.

# Introduction

What is ToPy?
ToPy is short for **t**opology **o**ptimisation (optimization) using **Py**thon.

ToPy can solve one of three types of topology optimisation problems: 
  1. minimum compliance (same as maximum stiffness),
  1. heat conduction or
  1. mechanism design (synthesis).

A problem is defined (created) by means of a simple text input file, called a TPD file, which stands for ToPy Problem
Definition.

ToPy solves the defined problem to obtain a 2D (or 3D, depending on the input file) checker-free black-white (in 2D) or solid-void (in 3D) solution. The result is
  1. an optimally stiff structure for minimum compliance (maximum stiffness) problems,
  1. an optimal distribution of two materials for heat conduction problems and
  1. an optimal distribution of material for efficient mobility.

The 2D results are PNG (or any format that Matplotlib can save) files and the 3D results are legacy VTK files.

There is a tutorial to solve a 2D compliance problem here: [Tutorials](Tutorials.md).

## Limitations
  * ToPy only works with regular square (for 2D)  and cubic (for 3D) meshes. Your mesh will therefore consist of perfectly square (for 2D) or perfectly cubic (for 3D) elements.
  * No GUI for defining problems
  * No CAD interface (although you can save the 3D files as STL files via ParaView, I'm sure there are other tricks too)

# Examples
See [Examples](Examples-of-ToPy-results.md) for some examples of the input (problem definition) files and their output after it was solved with ToPy.





