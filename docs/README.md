# Documentation

## Getting started
The main class of **ToPy** is 'Topology'. It defines the main constraints,
grid and parameters of optimization -- but you don't really have to bother
yourself with this if you just want to get some results.

### Defining a problem
There are two ways of defining a problem:
1. **TPD file**: You define the problem with keywords
(see [Help](docs/Help.md)) in a simple text file and solve via the command line. The text file must have the extension `.tpd`
2. **Config dictionary**: This is similar to the TPD file approach, however,
you define the problem directly in a Python file; it's very useful if you want to
experiment and don't want to keep making changes to a text file.
You can later save the Config keywords to a TPD file.

#### TPD (**T**oPy **P**roblem **D**efinition) file
There is a minimal set of parameters which is required for successful definition of a ToPy problem:
```
PROB_TYPE  : comp
PROB_NAME  : mbb_beam_minimal
ETA        : 0.5
DOF_PN     : 2
VOL_FRAC   : 0.5
FILT_RAD   : 1.5
P_FAC      : 3
ELEM_K     : Q4
NUM_ELEM_X : 60
NUM_ELEM_Y : 20
NUM_ELEM_Z : 0
NUM_ITER   : 10
FXTR_NODE_X: 1|21
FXTR_NODE_Y: 1281
LOAD_NODE_Y: 1
LOAD_VALU_Y: -1
```
You can read more about successful problem definition [here](./templates/).

When the TPD file is defined, then the rest is simple:

```python
from topy import Topology

topology = Topology()
topology.load_tpd_file('file.tpd')
```

#### Config dictionary
First you have to define a config dictionary (note the similarity with a TPD
file, especially the keywords):

```Python
config = {
     'DOF_PN': 2,
     'ELEM_K': 'Q4',
     'ETA': '0.5',
     'FILT_RAD': 1.5,
     'FXTR_NODE_X': range(1, 22),
     'FXTR_NODE_Y': 1281,
     'LOAD_NODE_Y': 1,
     'LOAD_VALU_Y': -1,
     'NUM_ELEM_X': 60,
     'NUM_ELEM_Y': 20,
     'NUM_ELEM_Z': 0,
     'NUM_ITER': 94,
     'PROB_NAME': 'beam_2d_reci',
     'PROB_TYPE': 'comp',
     'P_FAC': 3.0,
     'VOL_FRAC': 0.5
}
```
The requirements are the same as for the TPD file.

```Python
topology = Topology(config=config)
```
### Optimization (solving the problem)

You can use the command line solution:

```bash
$ python topy/scripts/optimise.py <filename>.tpd
```

Or you can use a Python script:

```Python
import topy

config = {...}
t = topy.Topology(config)
t.set_top_params()
topy.optimise(t)
```

The first time you run ToPy after a new install you'll see the
Element stiffness matrices creation in your terminal:

	It seems as though all or some of the element stiffness matrices
	do not exist. Creating them...
	This is usually only required once and may take a few minutes.
	SymPy is integrating: K for Q4bar...
	Created ToPy\topy\data\Q4bar.K (stiffness matrix).
	SymPy is integrating: K for Q4...
	Created ToPy\topy\data\Q4.K (stiffness matrix).
	SymPy is integrating: K for Q5B...
	Created ToPy\topy\data\Q5B.K (stiffness matrix).
	SymPy is integrating: K for Q4T...
	Created ToPy\topy\data\Q4T.K (stiffness matrix).
	SymPy is integrating: K for H8...
	Created ToPy\topy\data\H8.K (stiffness matrix).
	SymPy is integrating: K for H18B...
	Created ToPy\topy\data\H18B.K (stiffness matrix).
	SymPy is integrating: K for H8T...
	Created ToPy\topy\data\H8T.K (stiffness matrix).

You shouldn't see it again, even if ToPy is updated, since these
files shouldn't need to change. You can create the stiffness matrices without
solving a problem by simply running 'optimise.py' in the 'scripts' folder.


### Visualization (seeing the result)
Module `topy.visualization` allows one to save the output as a `.png` image for 2D problems or as a `.vtk` file for 3D.
The VTK files can be viewed with Mayavi or ParaView.
You can animate the PNG images with
the [convert](https://www.imagemagick.org/script/convert.php) tool.

```bash
convert -delay 35 *.png anim.gif
```

## Examples
See [Examples](Examples-of-ToPy-results.md) for some examples of the input (problem definition) files and their output after it was solved with ToPy.

## Tutorials
There is a tutorial to solve a 2D compliance problem here: [Tutorials](Tutorials.md).

## Original dissertation
The docs folder contains [William Hunter master's dissertation](thesis.pdf) since
ToPy was part of it. The maths behind topology optimization is explained in it.
There is also quite a bit about ToPy in it, but now a little outdated.
The persistent URL for the dissertation is: http://hdl.handle.net/10019.1/2648







