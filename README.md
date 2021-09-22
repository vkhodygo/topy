# ToPy
<div align="center">
	<img src="./imgsrc/topy_logo3d.png" width="400">
</div>

ToPy is a lightweight topology optimization framework for Python that can solve
compliance (stiffness), mechanism synthesis and heat conduction problems in 2D and 3D.

ToPy was originally created by [William Hunter](https://github.com/williamhunter/topy) for his Master's dissertation, in 2009. This fork was then created in 2020 so I could use it as a foundation to test an experimental stress-based sizing algorithm for topology optimization for my own Master's dissertation. I have both implemented the entire algorithm I developed and also upgraded the code for Python 3 and implemented whichever optimizations I could.

As Mr. Hunter actually did most of the work by fully implementing the topology optimization SIMP algorithm (and did a great job at it, I might add), I made sure to make the code backwards compatible. Not to mention that it also helps comparing both algorithms.

You can check Mr. Hunter's dissertation [here](http://hdl.handle.net/10019.1/2648).

Please refer to the [ToPy Wiki](https://github.com/TarcisioLOliveira/topy/wiki) for further information.

## Example of a ToPy result
An [example TPD file and solution/result](https://github.com/TarcisioLOliveira/topy/wiki/Examples-of-ToPy-results)

## Installation
Once you've downloaded the dependencies (see the [INSTALL](https://github.com/TarcisioLOliveira/topy/blob/master/INSTALL.md)
file) all you need to do is the following:

Download the latest **stable release** from here: https://github.com/TarcisioLOliveira/topy/releases/latest

Then do

```bash
cd topy/topy
python setup.py install
```

If you opted to use a virtual environment, make sure it is activated before and 
while you use ToPy.

## Getting started
The main classes of **ToPy** are 'TopologyGen' and 'TopologyTrad'.
They define the main constraints,
grid and parameters of optimization -- but you don't really have to bother
yourself with this if you just want to get some results.

### There are two ways of defining a problem
1. **TPD file**: You define the problem with keywords
(see [Help](https://github.com/TarcisioLOliveira/topy/wiki/Help)) in a simple text file and solve via the command line. The text file must have the extension `.tpd`
2. **Config dictionary**: This is similar to the TPD file approach, however,
you define the problem directly in a Python file. While it still works, I wouldn't
recommend it unless you are familiar with MPI.

### TPD (**T**oPy **P**roblem **D**efinition) file
There is a minimal set of parameters which is required for successful definition of a ToPy problem:
```
[ToPy Problem Definition File v2020]

TO_TYPE    : trad
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
To successfully use the new algorithm implemented, or take advantage of the stress
calculations for analysis, then the minimal set becomes:
```
[ToPy Problem Definition File v2020]
TO_TYPE    : gen         # or 'trad'

ELEM_E:      186         # GPa
ELEM_NU:     0.29
ELEM_TC:     51.9        # W/(m*K)
ELEM_L:      25          # mm (half-length of element)
THICKNESS:   50	         # mm
NORMAL_MAX : 100         # MPa
SHEAR_MAX  : 100         # MPa

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

Be aware that the new algorithm is only partially implemented for heat conduction,
and is not present for mechanism synthesis. Heat condunction also does not lead
to stress calculations (as it does not have a clear analogous concept).

You can read more about successful problem definition [here](https://github.com/TarcisioLOliveira/topy/tree/master/templates).

When the TPD file is defined, then the rest is simple. You can use the command line solution:

```bash
$ python topy/scripts/optimise.py <filename>.tpd
```

Or, for multiprocessing FEA (in this case, 2 processes): 

```bash
$ mpirun -n 2 python3 -m mpi4py optimise.py <filename>.tpd
```

### Config dictionary
First you have to define a config dictionary (note the similarity with a TPD
file, especially the keywords):

```Python
config = {
    'TO_TYPE': trad,
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
The requirements are the same as for the TPD file. Then do:

```Python
topology = TopologyTrad(config=config) # or TopologyGen, for 'TO_TYPE': gen
topology.set_top_params()
topy.optimise(t)
```

Just remember to make your code MPI aware.

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

For 3D problems with Paraview, for example, you can visualize them with:
```bash
paraview <filename>.vtk
```

## Tutorials
[Tutorials](https://github.com/TarcisioLOliveira/topy/wiki/Tutorials)

## How to cite ToPy
If you've used ToPy in your research work or find it useful in any way, please consider to cite:
```
@misc{TopyRepo
    author = {de Oliveira, Tarc\'isio Ladeia and Hunter, William and others},
    title = {ToPy - Topology optimization with Python},
    year = {2021},
    publisher = {GitHub},
    journal = {GitHub repository},
    howpublished = {\url{https://github.com/TarcisioLOliveira/topy}},
}
```
