# Dependencies
There always are... It shouldn't take more than a few minutes to download and install everything.

## All platforms
Install Python 3.x (32 or 64 bit) - make sure your 'bitness' is correct
when downloading the other packages listed further below.
Download Python from the official Python site.

## Linux
Install packages below via `pip` or by other means
(e.g., apt-get, yum, rpm,...). So, install he following:
1. NumPy+MKL
2. SciPy
3. matplotlib
4. SymPy
5. PyVTK

# Installing ToPy
CD into the 'topy' directory (where 'setup.py' is located) and
in a terminal, type:

	python setup.py install

or if you want to install locally, type:

	python setup.py install --user

If there aren't any errors, then ToPy is installed. Congratulations!

# Getting started
See [the tutorials](Tutorials.md)

# First run of ToPy
## Element stiffness matrices
The first time you run ToPy after a new install you'll see the
following in your terminal:

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

You won't (shouldn't) see it again, even if ToPy is updated, since these
files shouldn't need to change. You can create the stiffness matrices without
solving a problem by simply running 'optimise.py' in the 'scripts' folder.

---

William Hunter (2017)  
Francisco Sanchez Arroyo (2018)
