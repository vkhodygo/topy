![pytest](https://github.com/mlaradji/topy/workflows/pytest/badge.svg)

# ToPy
<div align="center">
	<img src="./imgsrc/topy_logo3d.png" width="400">
</div>

## What is ToPy?
ToPy is short for **t**opology **o**ptimization using **Py**thon. It is a lightweight framework for Python that can solve one of three types of topology optimisation problems: 
1. minimum compliance (same as maximum stiffness),
2. heat conduction or
3. mechanism design (synthesis).

## Example of a ToPy result
An [example TPD file and solution/result](https://github.com/williamhunter/ToPy/wiki/Examples-of-ToPy-results)

## Installation
**NOTE**: I've added a 0.4.1 release, which is older then the master branch, but will get you up and running with Python 2 and
Pysparse if you're willing to use the Anaconda Python distribution

Once you've downloaded the dependencies (see the [INSTALL](https://github.com/williamhunter/topy/blob/master/INSTALL.md)
file) all you need to do is the following:

Download the latest **stable release** from here: https://github.com/williamhunter/topy/releases/latest

ToPy solves the defined problem to obtain a 2D (or 3D, depending on the input file) checker-free black-white (in 2D) or solid-void (in 3D) solution. The result is
1. an optimally stiff structure for minimum compliance (maximum stiffness) problems,
2. an optimal distribution of two materials for heat conduction problems and
3. an optimal distribution of material for efficient mobility.

The 2D results are PNG (or any format that Matplotlib can save) files and the 3D results are legacy VTK files.


### ToPy and Python 3
ToPy is fairly old. I started working on it in 2005 and finished it around 2009, so that implies that the stable release only 
works with Python 2. You can however pull the latest "unstable" version, which should work with Python 3 (thanks to the
efforts of other people).

You will also need to install [Gmsh](http://gmsh.info).

Once you've downloaded the dependencies, all you need to do is the following:

```bash
$ git clone https://github.com/thebeachlab/topy.git
$ cd topy/topy
$ (sudo) python setup.py install
```
If there aren't any errors, then ToPy is installed. Congratulations!

## Documentation
Please refer to the [ToPy Docs](docs/README.md) for further information.

## Limitations
  * ToPy only works with regular square (for 2D)  and cubic (for 3D) meshes. Your mesh will therefore consist of perfectly square (for 2D) or perfectly cubic (for 3D) elements.
  * No GUI for defining problems
  * No CAD interface (although you can save the 3D files as STL files via ParaView)


### Visualization (seeing the result)
Module `topy.visualization` allows one to save the output as a `.png` image for 2D problems or as a `.vtk` file for 3D.
The VTK files can be viewed with Mayavi or ParaView.
You can animate the PNG images with
the [convert](https://www.imagemagick.org/script/convert.php) tool.

```bash
convert -delay 35 *.png anim.gif
```

<div align="left">
	<img src="./imgsrc/beam_2d_reci_gsf.gif" width=40%>
	<img src="./imgsrc/inverter_2d_eta03.gif" width=30%>
	<img src="./imgsrc/t-piece_2d_Q4_eta04_gsf.gif" width=20%>
</div>

## Tutorials
[Tutorials](https://github.com/williamhunter/topy/wiki/Tutorials)

## How to cite ToPy
If you've used ToPy in your research work or find it useful in any way, please consider to cite:
```
@misc{Hunter2007william,
  author = {Hunter, William and others},
  title = {ToPy - Topology optimization with Python},
  year = {2017},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/williamhunter/topy}},
  }
```
