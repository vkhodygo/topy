# ToPy Problem Definition keywords

This page explains what each ToPy keyword stands for as used in a TDP (ToPy Problem Definition) file. You should also look at the `template.tpd` file (in the source), which gives more information (also below).

## General keywords
| **KEYWORD**   | **Description**                                                                     |
| :------------ | :---------------------------------------------------------------------------------- |
| PROB\_NAME    | Problem name                                                                        |
| PROB\_TYPE    | Problem type                                                                        |
| NUM\_ELEM\_X  | Number of elements in the X-direction                                               |
| NUM\_ELEM\_Y  | Number of elements in the Y-direction                                               |
| NUM\_ELEM\_Z  | Number of elements in the Z-direction                                               |
| VOL\_FRAC     | Volume fraction                                                                     |
| FILT\_RAD     | Filter radius                                                                       |
| NUM\_ITER     | Number of iterations to run                                                         |
| CHG\_STOP     | Change stop value; controls the number of iterations                                |
| P\_FAC        | Start value of penalty factor (p)                                                   |
| P\_INCR       | Increment value of the penalty factor                                               |
| P\_CON        | Number of iterations to keep the penalty factor constant after each increment       |
| P\_MAX        | Max value that the penalty factor is allowed to reach                               |
| Q\_FAC        | Start value of extra penalty factor (q) for grey-scale filtering (GSF)              |
| Q\_INCR       | Analoguous to P\_INCR                                                               |
| Q\_CON        | Analoguous to P\_CON                                                                |
| Q\_MAX        | Analoguous to P\_MAX                                                                |
| ETA           | Damping factor                                                                      |
| APPROX        | Type of approximation, reciprocal, exponential or diagonal quadratic                |
| ELEM\_K       | Type of finite element to use; really specifies the element stiffness matrix to use |
| FXTR\_NODE\_X | Node number(s) at which you want to fix the translation in the X direction          |
| FXTR\_NODE\_Y | Node number(s) at which you want to fix the translation in the Y direction          |
| FXTR\_NODE\_Z | Node number(s) at which you want to fix the translation in the Y direction          |
| LOAD\_NODE\_X | Node number(s) at which you require load(s) in the X direction                      |
| LOAD\_NODE\_Y | Node number(s) at which you require load(s) in the Y direction                      |
| LOAD\_NODE\_Z | Node number(s) at which you require load(s) in the Z direction                      |
| LOAD\_VALU\_X | Load value(s) in the X direction                                                    |
| LOAD\_VALU\_Y | Load value(s) in the Y direction                                                    |
| LOAD\_VALU\_Z | Load value(s) in the Z direction                                                    |
| PASV\_ELEM    | Element number(s) at which you require passive (void) elements                      |
| ACTV\_ELEM    | Element number(s) at which you require active (solid) elements                      |

### Mechanism synthesis
| **KEYWORD**        | **Description**                                                   |
| :----------------- | :---------------------------------------------------------------- |
| LOAD\_NODE\_X\_OUT | Node number(s) at which you require the output in the X direction |
| LOAD\_VALU\_X\_OUT | Value of output at specified OUT node in X direction              |
| LOAD\_NODE\_Y\_OUT | Node number(s) at which you require the output in the Y direction |
| LOAD\_VALU\_Y\_OUT | Value of output at specified OUT node in Y direction              |
| LOAD\_NODE\_Z\_OUT | Node number(s) at which you require the output in the Z direction |
| LOAD\_VALU\_Z\_OUT | Value of output at specified OUT node in Z direction              |

***


# TPD file format
The file *must* start with the identifier `[ToPy Input File v2007]` and then
 a blank line.

Comments may be placed after the the hash \# character, in-line comments are also supported.

Order of keywords (parameters) may be random, but keywords need to be first
 in sequence, e.g:
* KEYWORD1: some value
* KEYWORD2: another value
* KEYWORD3: yet another value
* KEYWORD4: etc.

There is no restriction on the use of whites pace, but TABs will create problems, thus, make sure your editor replaces TABs with SPACEs (the ToPy parser will warn you if it finds TABs anyway).


## Problem types
  Can be one of the following
* 'comp' = minimum compliance problem,
* 'heat' = heat conduction problem or
* 'mech' = mechanism design (syntehesis) problem.

The keyword-value case doesn't matter, i.e., comp = CoMp = COMP.
    
    PROB_TYPE: comp # Solve a minimum compliance problem.

Problem name:

    PROB_NAME: really_cool_problem # Output files will have this name.


## Problem parameters (keywords)

    VOL_FRAC: 0.5 # The volume fraction to be used.
    FILT_RAD: 1.5 # Filter radius.

Use one of the following:

    NUM_ITER: 100 # Number of iterations to run.
    CHG_STOP: 0.01 # Change stop value, checks the change in obj. function value.

    P_FAC : 3      # Start value of penalty factor (p).
    P_MAX : 3.5    # Max value of P_FAC.
    P_INCR: 0.02   # Increment value of P_FAC.
    P_CON : 25     # Number of iterations to keep P_FAC constant after increment.
    Q_FAC : 1      # Start value of extra penalty factor (q) for GSF.
    Q_MAX : 5      # Analoguous to P_MAX.
    Q_INCR: 0.08   # Analoguous to P_INCR.
    Q_CON : 20     # Analoguous to P_CON.

    ETA   : 0.5     # 0.001, use 0.5 for reciprocal approximation.
    ETA   : exp     # Use exponential approximation, eta is 'auto-tuned'
    APPROX: dquad   # Use diagonal quadratic approximation, ETA must be specified.

## Finite Element Types
Nodes are numbered as follows, although this is not important to the user.
```
 2D: Y             3D: Y
     |                 |
   4-|-3             4-|-3
   | +-|---X        /| +-|---X
   1---2           / 1/--2
                  8--/7 /
                  | / |/
                  5/--6
                  /
                 Z

ELEM_K: Q4 # Other 2D: Q5B, Q4a5B, Q4T.   3D: H8, H18B, H8T.
```

## Discretisation of the design domain
```
 2D: Y             3D: Y
     |                 |
     +---X             +---X
                      /
                     Z

 1---5---9
 | 1 | 5 |
 2---6---10
 | 2 | 6 |
 3---7---11
 | 3 | 7 |
 4---8---12
```

Numbering of nodes and elements is from top to bottom, column-wise, starting
 at one (1). For 3D, the X-Y plane is numbered first, then in the Z-direction.

    NUM_ELEM_X: 60 # Number (quantity) of elements in the X-direction.
    NUM_ELEM_Y: 20 # Number of elements in the (negative) Y-direction.

Set the following keyword to 0 if not necessary for your problem, i.e., 2D:
    NUM_ELEM_Z: 10

## Translational constraints
Node number(s) and/or 'start|stop|step' notation may be used for multiple
 nodes, ";" (semi-colon) may be used to separate ranges.
NOTE: Do not end a line with a semi-colon nor use commas anywhere.

FXTR_NODE_X = FiX TRanslation of NODE in the X direction

    FXTR_NODE_X: 1|21               # Node 1 to 21, step size 1 is implied.
    FXTR_NODE_Y: 1281               # Lower right corner for 60x20 problem
    FXTR_NODE_Z: 1; 4|13|3; 18|22|2 # Nodes 1, 4, 7, 10, 13, 18, 20, 22.

## Loads
 Node number(s) and/or 'start|stop|step' notation may be used for multiple
 nodes, ";" may be used to separate ranges. NOTE: Do not end a line with a ";"
 A load value (VALU) *must* be specified for the nodes you choose. Use + or -
 to set the direction of the load.
 Set the *node number(s)* that's loaded, *not* the degrees of freedom, that's
 taken care of by ToPy. Assign values for the corresponding load size(s).
 Also note the use of the "@" below, which is rather convenient for the user.
    
    LOAD_NODE_X: 1; 4; 9 # Load nodes 1, 4 and 9 in X direction.
    LOAD_NODE_Y: 1       # Upper left corner -- node number = 1 (always), Y direction.
    LOAD_NODE_Z: 20|32|3 # Load node 20 to 32 in steps of 3, in Z direction.
    LOAD_VALU_X: 0.75    # Simply omit a line if not necessary for your problem.
    LOAD_VALU_Y: -1      # Value of the load = 1, direction negative Y.
    LOAD_VALU_Z: 1@10    # Value of the load = 1 at 10 nodes in Z direction.


## Passive (void) and active (solid) elements
List the *element* numbers you want to affect.
    
    PASV_ELEM: 10|19; 30|39; 50|59; 70|79; 90|99 # No elements will appear at these locations (void)
    ACTV_ELEM: 1|1181|20; 1181|1200 # These elements will remain (solid)


## Mechanism design (synthesis) specific

    LOAD_NODE_X_OUT: 841 # Node number(s) at which you require the output
    LOAD_VALU_X_OUT: -1 # Value of output at specified OUT node.


## Heat conduction specific
NOTE: For heat conduction problems, only use *_X keywords, i.e., no Y or Z
 dimensions, since heat problems are one-dimensional i.t.o. degrees of freedom
 (temperature is a scalar value).
However, you still have to specify the following for 2D problems:

    NUM_ELEM_Z : 0

And don't forget this, for example:

    DOF_PN  : 1
    ELEM_K  : Q4T

That's it, easy peasy :-)