# Dependencies
This version features more dependencies than the original one, but most are
relatively easy do get (meaning, by using pip).

The main dependency is the sparse matrix direct solver [MUMPS](http://mumps.enseeiht.fr/index.php?page=home).
For this package, you have to install MUMPS along with MPI (or OpenMPI). On Linux
and Mac you can install it more easily using [Homebrew/Linuxbrew](https://brew.sh/).
Install it with:

```bash
brew tap brewsci/num
brew install brewsci-mumps
```

Make sure you have both a C and a FORTRAN compiler installed. GCC has both by
default. 

On Windows, you may have to compile MUMPS yourself. [This repo](https://github.com/scivision/mumps)
maintains a fork which allows you to build MUMPS using CMake, and may be the
better option. Just remember you need to build it with MPI support.

Other dependencies are Python packages, you may get them with `pip`:

```bash
pip install typing pathlib matplotlib sympy numpy pyvtk scipy cython mpi4py pymumps numexpr
```

Be aware that you may need to point to MPI's and MUMPS's include and binary files
to be able to install `mpi4py`and `pymumps`. If you installed them using Linuxbrew,
they will be at:

    /home/linuxbrew/.linuxbrew/opt/brewsci-mumps/lib
    /home/linuxbrew/.linuxbrew/opt/brewsci-mumps/include
    /home/linuxbrew/.linuxbrew/opt/openmpi/lib
    /home/linuxbrew/.linuxbrew/opt/openmpi/include
    /home/linuxbrew/.linuxbrew/opt/openmpi/bin

Add the last one to your `PATH` with:

```bash
export PATH="/home/linuxbrew/.linuxbrew/opt/openmpi/bin:$PATH"
```

And then run, in this order:

```bash
pip install mpi4py --global-option="build_ext" \
  --global-option="-I/home/linuxbrew/.linuxbrew/opt/openmpi/include" \
  --global-option="-L/home/linuxbrew/.linuxbrew/opt/openmpi/lib"
pip install pymumps --global-option="build_ext" \
  --global-option="-I/home/linuxbrew/.linuxbrew/opt/brewsci-mumps/include" \
  --global-option="-L/home/linuxbrew/.linuxbrew/opt/brewsci-mumps/lib"
```

It is a good idea to keep OpenMPI's `bin`folder in your path to run ToPy, especially
because you need `mpiexec` to run MUMPS in parallel. To do that, I recommend using
a virtual environment and include them in its `activate` file. For example:

```bash
python3 -m venv ~/venv-topy
gedit ~/venv-topy/bin/activate
```

And then add:

    PATH="/home/linuxbrew/.linuxbrew/opt/openmpi/bin:$PATH"

After the function `deactivate ()`. You can then enter the virtual environment
with `source ~/venv-topy/bin/activate` and leave it with `deactivate`.

# Installing ToPy
Download the stable topy release or clone.

CD into the 'topy' directory (where 'setup.py' is located) and
in a terminal ('cmd' on Windows), type:

```bash
python setup.py install
```

The command above is the one to use if you are using a virtual environment too.
If you want to install locally, type:

```bash
python setup.py install --user
```

You may require Administrator rights on Windows, depending on your setup.

If there aren't any errors, then ToPy is installed. Congratulations!

# Getting started
See https://github.com/williamhunter/topy/wiki/Tutorials

# First run of ToPy
The original ToPy created all element matrices on the first run and stored them
for future use. This version introduces the possibility to use material properties
such as Young's modulus and thermal conductivity, so now the element required by
the TPD file will be generated with each execution.
