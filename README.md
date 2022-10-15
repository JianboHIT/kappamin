# kappamin (under-developed)
A python3 code for calculations of the minimum limit to thermal conductivity

## Features

- Models of the minimum limit to thermal conductivity under Cahill assumption
  - Debye model
  - BvK (Born–von Karman) model
  - Pei model
- Temperature-dependence
  - Finite temperature
  - Ideal infinite temperature
- Relative
  - Heat Capacity
  - Minimum mean-free-path
  - Minimum average phonon lifetime
- Running mode
  - Command line mode based on a simple configuration file (for the routine analysis)
  - Based on prepared scripts (for general researchers without programming skills)
  - Use as a python module (for expert usage)

# Getting Started

In order to install the module, you need to install the dependencies first:

- python3
- numpy
- scipy

Download the package from GitHub website or using `git`,

```
$ git clone https://github.com/JianboHIT/kappamin.git
```

then run `setup.py` to install.

```
$ cd kappamin
$ python3 setup.py install
```

Python module `kappamin` will be installed if no error. 
In order to invoke `kappamin` module, you need to prepare a configuration file 
(see [Example_Debye.txt](Example_Debye.txt)
and [Example_BvK.txt](Example_Debye.txt)
in the source package).

```
$ python -m kappamin [KAPPAMIN.txt]
```

Here `KAPPAMIN.txt` indicates the filename of configuration file. 
It is worth mentioning that the filename is optional. 
If the filename is not given, the program will read the file named as `KAPPAMIN.txt` if it existed.

Alternately, a more convenient way to implement calculation is by an executable script
(see [ExceuteScript.py](ExecuteScript.py)), then run it by python3. 
On Linux or Windows Terminal:

```
$ python ExceuteScript.py
```

On Windows, if it has been configured that the default program to open .py file is python3,
you just need to move ExceuteScript.py to the directory at where the configuration file is located
and double-click it to run.

## Feedback and report bugs

See [GitHub Issue page](../../issues).

## Change log

- 2022.10.06 v0.0.1 Initial package version

