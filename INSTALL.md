# Prerequisites

1. Start by installing Python 2.7 (32-bit works fine), download from the official site. If you're using Linux, you most probably have it already.
1. **Windows only:** Make sure Python is in your Path Environment Variable.
	1. If you don't know how to add it, please search the web.
	1. Check if Python works by typing it into a *cmd* shell.
	1. Check if 'pip' works, also by typing it into a *cmd* shell. If it doesn't work, add Python27\Scripts to the Environment Variables too.
	1. When using 'pip', run 'cmd' as Administrator.
	1. Install NumPy+MKL for Python 2.7 from here: http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy
	1. Install PySparse, also from http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy
	1. Install PyVTK, also from http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy
1. **Linux**:
	1. On Ubuntu 16.04, run the following commands in a terminal: (they should normally work on any debian-based OS)
		1. `sudo apt-get install -y python-pip python-dev python-tk libpython-dev libevent-dev libsuperlu-dev libblas-dev liblapack-dev libatlas3-base libatlas-dev`
		2. `pip install --upgrade pip`
		3. `sudo pip install matplotlib setuptools SymPy pysparse pyvtk`
		4. Then download manually the files you want or download the full project (run in a terminal the following: `git clone https://github.com/williamhunter/topy`)
		5. According to what you are waiting from the software, if you want to do 3D topology optimizations, then you should also install Paraview: `sudo apt-get install -y paraview`
	1. On other Linux OS, you should install above packages (NumPy+MKL, Pysparse, PyVTK) via `pip` or by other means (e.g., apt-get, yum, rpm, etc). Then:
		1. Install matplotlib via `pip`
		1. Install SymPy via `pip`

Installing matplotlib and SymPy via other 'official' channels should also work fine (in that ToPy should still work).

If everything installed correctly, you're set.

# Installing ToPy
In a terminal ('cmd' window on Windows), type:

	python setup.py install

or if you want to install locally, type:

	python setup.py install --user

You may require Administrator rights on Windows, depending on your setup.
	
## Creating the element stiffness matrices
The first time you run ToPy after a fresh install you'll see the following in your terminal:

	It seems as though all or some of the element stiffness matrices
	do not exist. Creating them...
	This is usually only required once and may take a few minutes.
	SymPy is integrating: K for Q4bar...
	Created C:\Users\William\Programming\ToPy\topy\data\Q4bar.K (stiffness matrix).
	SymPy is integrating: K for Q4...
	Created C:\Users\William\Programming\ToPy\topy\data\Q4.K (stiffness matrix).
	SymPy is integrating: K for Q5B...
	Created C:\Users\William\Programming\ToPy\topy\data\Q5B.K (stiffness matrix).
	SymPy is integrating: K for Q4T...
	Created C:\Users\William\Programming\ToPy\topy\data\Q4T.K (stiffness matrix).
	SymPy is integrating: K for H8...
	Created C:\Users\William\Programming\ToPy\topy\data\H8.K (stiffness matrix).
	SymPy is integrating: K for H18B...
	Created C:\Users\William\Programming\ToPy\topy\data\H18B.K (stiffness matrix).
	SymPy is integrating: K for H8T...
	Created C:\Users\William\Programming\ToPy\topy\data\H8T.K (stiffness matrix).
	
You won't (shouldn't) see it again, even if ToPy is updated, since these
files shouldn't need to change.
